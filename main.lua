-- craquelure viewer for LOVE 11.4+ (love2d.org)
-- runs the phase-field film fracture model step by step and shows the
-- crack network forming while the shrinkage mismatch (load) grows

local craq = require("craquelure")
local presets = require("data_presets")

local state = {
  preset = 1,
  size_index = 2,
  sizes = { 128, 192, 256, 384, 512 },
  seed = 1,
  view = 1,
  views = { "craquelure", "damage", "stress", "thickness" },
  scale = 2,           -- display image has n*scale pixels per side
  export_scale = 4,
  precise = false,     -- staggered tolerance 1e-3 instead of 5e-3
  width_scale = 0.5,   -- crack line width per unit of crack opening
  tile = false,
  paused = false,
  budget = 0.030,      -- solver time per frame, seconds
  solver = nil,
  co = nil,
  imagedata = nil,
  image = nil,
  quad = nil,
  dirty = true,
  last_refresh = 0,
  last_step = -1,
  message = "",
  message_time = 0,
  finished_at = nil,
  started_at = 0,
}

-- coroutine body: runs the whole quasi-static loading
local function run_body(S)
  craq.run(S)
  return "done"
end

local function copy_params(src)
  local t = {}
  for k, v in pairs(src) do
    t[k] = v
  end
  return t
end

local function start()
  if state.solver then
    state.solver.stop = true
  end
  state.solver = nil
  state.co = nil
  state.imagedata = nil
  state.image = nil
  collectgarbage()
  collectgarbage()
  local pr = presets[state.preset]
  local p = copy_params(pr.params)
  p.n = state.sizes[state.size_index]
  p.seed = state.seed
  p.tol_alpha = state.precise and 1e-3 or 5e-3
  -- grid size changes the size of the film patch, not the cell size
  local S = craq.new(p)
  S.yield_fn = coroutine.yield
  state.solver = S
  state.co = coroutine.create(run_body)
  local W = S.N * state.scale
  state.imagedata = love.image.newImageData(W, W)
  state.image = love.graphics.newImage(state.imagedata)
  state.image:setFilter("linear", "linear")
  state.image:setWrap("repeat", "repeat")
  state.quad = love.graphics.newQuad(0, 0, 2 * W, 2 * W, W, W)
  state.dirty = true
  state.last_step = -1
  state.isl_count, state.isl_area = nil, nil
  state.finished_at = nil
  state.started_at = love.timer.getTime()
end

local function refresh_image()
  local S = state.solver
  if not S then
    return
  end
  craq.render(S, state.imagedata:getFFIPointer(), state.views[state.view], state.scale,
    state.width_scale)
  state.image:replacePixels(state.imagedata)
  state.dirty = false
  state.last_refresh = love.timer.getTime()
end

local function flash(text)
  state.message = text
  state.message_time = love.timer.getTime()
end

local function export_pngs()
  local S = state.solver
  local W = S.N * state.export_scale
  local names = {}
  local modes = { "craquelure", "damage" }
  for _, mode in ipairs(modes) do
    local img = love.image.newImageData(W, W)
    craq.render(S, img:getFFIPointer(), mode, state.export_scale, state.width_scale)
    local name = string.format("craquelure_p%d_seed%d_step%d_%s_%dpx.png",
      state.preset, state.seed, S.step, mode, W)
    img:encode("png", name)
    names[#names + 1] = name
  end
  flash("exported " .. table.concat(names, ", ") .. " to " .. love.filesystem.getSaveDirectory())
end

function love.load()
  love.graphics.setBackgroundColor(0.11, 0.10, 0.09)
  start()
end

function love.update(dt)
  local co = state.co
  local S = state.solver
  if co and not state.paused and coroutine.status(co) ~= "dead" then
    local t_end = love.timer.getTime() + state.budget
    local running = true
    while running do
      local ok, err = coroutine.resume(co, S)
      if not ok then
        error("solver error: " .. tostring(err))
      end
      if coroutine.status(co) == "dead" then
        running = false
        state.finished_at = love.timer.getTime()
        state.dirty = true
      elseif love.timer.getTime() >= t_end then
        running = false
      end
    end
    if S.step ~= state.last_step or love.timer.getTime() - state.last_refresh > 1.0 then
      if S.step ~= state.last_step then
        state.isl_count, state.isl_area = craq.islands(S)
      end
      state.last_step = S.step
      state.dirty = true
    end
  end
  if state.dirty then
    refresh_image()
  end
end

local function draw_panel(x, y, w, h)
  love.graphics.setColor(0, 0, 0, 0.62)
  love.graphics.rectangle("fill", x, y, w, h, 6, 6)
  love.graphics.setColor(1, 1, 1, 1)
end

function love.draw()
  local Wn, Hn = love.graphics.getDimensions()
  local S = state.solver
  if not S or not state.image then
    return
  end
  local W = S.N * state.scale
  local side = math.min(Wn - 380, Hn - 20)
  if side < 200 then
    side = math.min(Wn, Hn) - 20
  end
  local ox, oy = 10, 10
  love.graphics.setColor(1, 1, 1, 1)
  if state.tile then
    love.graphics.draw(state.image, state.quad, ox, oy, 0, side / (2 * W), side / (2 * W))
  else
    love.graphics.draw(state.image, ox, oy, 0, side / W, side / W)
  end
  -- info panel
  local px = ox + side + 10
  local pw = Wn - px - 10
  if pw < 260 then
    px = ox
    pw = 360
  end
  draw_panel(px, 10, pw, Hn - 20)
  local pr = presets[state.preset]
  local lines = {}
  lines[#lines + 1] = string.format("[%d] %s", state.preset, pr.name)
  lines[#lines + 1] = ""
  for _, l in ipairs(pr.info) do
    lines[#lines + 1] = l
  end
  lines[#lines + 1] = ""
  lines[#lines + 1] = string.format("grid %d x %d   seed %d   %s", S.N, S.N, state.seed,
    state.precise and "precise" or "fast")
  lines[#lines + 1] = string.format("ell %.2f   le %.2f   le/ell %.1f", S.ell, S.le, S.le / S.ell)
  if S.path then
    lines[#lines + 1] = string.format("load path: state %d of %d", math.max(S.path_index, 1), #S.path)
  else
    lines[#lines + 1] = string.format("load p/p_ref %.3f  (to %.2f)", S.p / S.pref, S.load_end)
  end
  lines[#lines + 1] = string.format("mismatch %.2f / %.2f at %.0f deg", S.m1, S.m2, S.mangle)
  lines[#lines + 1] = "  (principal values, 1 = onset under equal shrinkage)"
  lines[#lines + 1] = string.format("load step %d   iterations in step %d", S.step, S.iter_in_step or 0)
  if state.isl_count and state.isl_count > 1 then
    local size = math.sqrt(state.isl_area)
    lines[#lines + 1] = string.format("islands %d, mean size %.0f cells = %.1f le",
      state.isl_count, size, size / S.le)
  end
  local status
  if state.finished_at then
    status = string.format("done in %.0f s", state.finished_at - state.started_at)
  elseif state.paused then
    status = "paused"
  else
    status = string.format("running %.0f s", love.timer.getTime() - state.started_at)
  end
  lines[#lines + 1] = status
  lines[#lines + 1] = "view: " .. state.views[state.view]
  lines[#lines + 1] = ""
  lines[#lines + 1] = "1-" .. #presets .. "  preset"
  lines[#lines + 1] = "space  pause"
  lines[#lines + 1] = "r  restart with new seed"
  lines[#lines + 1] = "v  view: craquelure / damage / stress / thickness"
  lines[#lines + 1] = "t  2x2 tiles (the pattern is periodic)"
  lines[#lines + 1] = "+ -  grid size (larger or smaller film patch)"
  lines[#lines + 1] = "a  precise / fast convergence"
  lines[#lines + 1] = ", .  crack line width per opening (" .. string.format("%.2f", state.width_scale) .. ")"
  lines[#lines + 1] = "s  save view as png"
  lines[#lines + 1] = "e  export craquelure and damage maps x" .. state.export_scale
  if love.timer.getTime() - state.message_time < 6 then
    lines[#lines + 1] = ""
    lines[#lines + 1] = state.message
  end
  love.graphics.printf(table.concat(lines, "\n"), px + 10, 20, pw - 20)
  -- load progress bar
  local frac = math.min(math.max(S.progress or 0, 0), 1)
  love.graphics.setColor(1, 1, 1, 0.25)
  love.graphics.rectangle("fill", px + 10, Hn - 34, pw - 20, 8)
  love.graphics.setColor(0.95, 0.75, 0.4, 1)
  love.graphics.rectangle("fill", px + 10, Hn - 34, (pw - 20) * frac, 8)
  love.graphics.setColor(1, 1, 1, 1)
end

function love.keypressed(key)
  if key == "escape" then
    love.event.quit()
  elseif key == "space" then
    state.paused = not state.paused
  elseif key == "r" then
    state.seed = state.seed + 1
    start()
  elseif key == "v" then
    state.view = state.view % #state.views + 1
    state.dirty = true
  elseif key == "t" then
    state.tile = not state.tile
  elseif key == "a" then
    state.precise = not state.precise
    start()
  elseif key == "=" or key == "+" or key == "kp+" then
    if state.size_index < #state.sizes then
      state.size_index = state.size_index + 1
      start()
    end
  elseif key == "-" or key == "kp-" then
    if state.size_index > 1 then
      state.size_index = state.size_index - 1
      start()
    end
  elseif key == "," then
    state.width_scale = math.max(state.width_scale - 0.25, 0)
    state.dirty = true
  elseif key == "." then
    state.width_scale = math.min(state.width_scale + 0.25, 3)
    state.dirty = true
  elseif key == "s" then
    local name = string.format("craquelure_p%d_%s_seed%d_step%d.png", state.preset,
      state.views[state.view], state.seed, state.solver.step)
    state.imagedata:encode("png", name)
    flash("saved " .. love.filesystem.getSaveDirectory() .. "/" .. name)
  elseif key == "e" then
    export_pngs()
  else
    local k = tonumber(key)
    if k and presets[k] then
      state.preset = k
      start()
    end
  end
end
