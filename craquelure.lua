-- craquelure.lua
-- physically based craquelure: a brittle film (paint, varnish, glaze, gesso)
-- bonded to a substrate cracks because it wants to shrink more than the
-- substrate lets it (drying, aging, cooling, humidity).
--
-- model, per unit film area, lengths in grid cells, young modulus E = 1:
--
--   E(u, a) = int h * [ g(a) * psi_plus(e) + psi_minus(e) ]   film elastic energy
--           + int k/2 * |u|^2                                 shear-lag foundation
--           + int h * gc * 3/8 * ( a/l + l * |grad a|^2 )     fracture energy (AT1)
--
--   u        in-plane film displacement relative to the substrate
--   e        elastic strain = sym(grad u) + m
--   m        mismatch strain: how much more the film wants to shrink than
--            the substrate lets it, a tensor (principal values and their
--            direction). either m = p * S with a growing load p and a fixed
--            shape S (s1, s2, angle), or any load path m(t), for example
--            drying shrinkage followed by humid and dry seasons
--   a        damage, 0 = intact, 1 = broken (smeared crack of width ~4 l)
--   g(a)     (1 - eta) * (1 - a)^2 + eta
--   psi+/-   plane stress energy split into tension and compression parts
--            (spectral split), so the film never cracks in compression
--   h        local film thickness factor (mean 1)
--   k        ebar / (le^2 * h^kexp), le = shear-lag length of the unit film;
--            with kexp = 1 this reproduces the exact energy release rate of a
--            channel crack in a film on a substrate, G = Z*sigma^2*h/ebar,
--            when le = Z*h (Z = 1.976 for equal film and substrate moduli)
--   l        AT1 length, a material constant: it sets the tensile strength
--            sigma_c = sqrt(3*gc*E/(8*l)) (uniaxial)
--
-- open cracks: where the damage reaches a_open the film is broken through.
-- such points are set to a = 1 and no longer bonded to the substrate (there
-- is no film left there to hold), so a crack band does not widen while
-- the crack opens.
--
-- evolution: quasi-static, the mismatch changes in small steps; at every
-- step the energy is minimized under irreversibility a >= a_previous_step by
-- alternate minimization: u by a newton step solved with jacobi pcg,
-- a by projected gauss-seidel. while cracks run, only the region around the
-- changing damage is re-solved (the film reacts within a few le), and every
-- load step ends with a global check. the grid is periodic, so results tile.

local ffi = require("ffi")
local bit = require("bit")

pcall(ffi.cdef, "void *malloc(size_t size); void free(void *ptr);")

local sqrt, abs, floor, ceil = math.sqrt, math.abs, math.floor, math.ceil
local max, min = math.max, math.min
local exp, log, cos, sin, pi = math.exp, math.log, math.cos, math.sin, math.pi

local G1 = 0.5 - 0.5 / sqrt(3)
local G2 = 0.5 + 0.5 / sqrt(3)
local Q11 = G1 * G1
local Q12 = G1 * G2
local Q22 = G2 * G2
local W = 0.25

local M = {}

M.defaults = {
  n = 256,            -- grid nodes per side, periodic
  ell = 2.0,          -- AT1 material length in cells (sets the strength)
  le = 8.0,           -- shear-lag length of the unit-thickness film, cells
  nu = 0.3,           -- poisson ratio of the film
  gc = 1.0,           -- fracture energy, units of E * cell
  eta = 1e-6,         -- residual stiffness of broken material
  kexp = 1.0,         -- foundation k ~ 1/h^kexp: 1 = film on substrate, 0 = soft interlayer
  split = true,       -- tension / compression split
  s1 = 1.0,           -- principal shrinkage along angle
  s2 = 1.0,           -- principal shrinkage across angle (negative = compression)
  angle = 0.0,        -- direction of s1, degrees
  gc_std = 0.1,       -- relative scatter of fracture energy (log-normal)
  gc_corr = 2.0,      -- correlation length of that scatter, cells
  h_std = 0.0,        -- relative scatter of thickness
  h_corr = 32.0,      -- correlation length of thickness, cells
  flaws = 0,          -- number of pre-existing micro-cracks
  flaw_len = 6.0,     -- their length, cells
  load_start = 0.9,   -- first load, in units of p_ref (homogeneous damage onset)
  load_end = 2.0,     -- last load, same units
  load_step = 0.01,   -- load increment, same units
  path = false,       -- optional load path instead of the load range, see M.run
  tol_alpha = 5e-3,   -- staggered convergence: max damage change per iteration
  max_iter = 20000,   -- staggered iterations per load step
  pgs_sweeps = 4,     -- symmetric gauss-seidel sweeps per damage update
  omega = 1.5,        -- over-relaxation of the damage update (1 = plain)
  a_open = 0.9,       -- damage above this counts as an open crack (a = 1):
                      -- the stiffness left there is below 1 percent
  cg_rtol = 1e-2,     -- relative pcg tolerance per newton step
  cg_maxit = 400,
  tile = 8,           -- tile size for the local solves, cells
  local_radius = 2.5, -- radius of a local solve, in shear-lag lengths
  local_frac = 0.5,   -- go global when more than this part of the film is active
  seed = 1,
}

-- memory -------------------------------------------------------------------

local function new_array(n, ctype)
  local size = ffi.sizeof(ctype)
  local raw = ffi.C.malloc(n * size)
  if raw == nil then
    error("craquelure: out of memory")
  end
  local a = ffi.cast(ctype .. "*", raw)
  ffi.fill(a, n * size, 0)
  return ffi.gc(a, ffi.C.free)
end

-- random numbers (xorshift32, same sequence everywhere) ----------------------

local function rng_new(seed)
  local s = bit.tobit(seed * 2654435761 + 12345)
  if s == 0 then
    s = 0x2545F491
  end
  return { s = s }
end

local function rng_uniform(r)
  local x = r.s
  x = bit.bxor(x, bit.lshift(x, 13))
  x = bit.bxor(x, bit.rshift(x, 17))
  x = bit.bxor(x, bit.lshift(x, 5))
  r.s = x
  return ((x % 4294967296) + 0.5) / 4294967296
end

local function rng_normal(r)
  local u1 = rng_uniform(r)
  local u2 = rng_uniform(r)
  return sqrt(-2 * log(u1)) * cos(2 * pi * u2)
end

-- correlated gaussian field, periodic, zero mean, unit variance
local function gaussian_field(N, corr, r, out, tmp)
  local n = N * N
  for i = 0, n - 1 do
    out[i] = rng_normal(r)
  end
  if corr > 0.3 then
    local rad = ceil(3 * corr)
    if rad > floor(N / 2) - 1 then
      rad = floor(N / 2) - 1
    end
    local kw = new_array(2 * rad + 1, "double")
    local ksum = 0
    for d = -rad, rad do
      local v = exp(-0.5 * d * d / (corr * corr))
      kw[d + rad] = v
      ksum = ksum + v
    end
    for d = 0, 2 * rad do
      kw[d] = kw[d] / ksum
    end
    for j = 0, N - 1 do
      local row = j * N
      for i = 0, N - 1 do
        local s = 0
        for d = -rad, rad do
          s = s + kw[d + rad] * out[row + (i + d) % N]
        end
        tmp[row + i] = s
      end
    end
    for j = 0, N - 1 do
      for i = 0, N - 1 do
        local s = 0
        for d = -rad, rad do
          s = s + kw[d + rad] * tmp[((j + d) % N) * N + i]
        end
        out[j * N + i] = s
      end
    end
  end
  local mean = 0
  for i = 0, n - 1 do
    mean = mean + out[i]
  end
  mean = mean / n
  local var = 0
  for i = 0, n - 1 do
    local d = out[i] - mean
    var = var + d * d
  end
  local sd = sqrt(var / n)
  if sd < 1e-300 then
    sd = 1
  end
  for i = 0, n - 1 do
    out[i] = (out[i] - mean) / sd
  end
end

-- constitutive law at one gauss point ----------------------------------------
-- elastic strain (exx, eyy, gxy engineering shear), pw = degraded weight,
-- mw = intact weight. stores W * tangent (c11 c22 c33 c12 c13 c23) at cg[o],
-- returns W * stress and the unweighted energies psi+ and psi-.
-- branch free on purpose (min / max only): branches inside these hot
-- functions make luajit fall back to the interpreter

local function gp_split(cg, o, exx, eyy, gxy, pw, mw, lam, mu, split)
  local sxx, syy, sxy, psip, psim
  if split then
    local tr = exx + eyy
    local hd = 0.5 * (exx - eyy)
    local exy = 0.5 * gxy
    local r = sqrt(hd * hd + exy * exy)
    local m = 0.5 * tr
    local e1 = m + r
    local e2 = m - r
    local rr = max(r, 1e-300)
    local c2 = hd / rr
    local s2 = exy / rr
    c2 = c2 + (1 - c2 * c2 - s2 * s2)
    local tp = max(tr, 0)
    local tn = min(tr, 0)
    local e1p = max(e1, 0)
    local e2p = max(e2, 0)
    local e1n = min(e1, 0)
    local e2n = min(e2, 0)
    psip = 0.5 * lam * tp * tp + mu * (e1p * e1p + e2p * e2p)
    psim = 0.5 * lam * tn * tn + mu * (e1n * e1n + e2n * e2n)
    local tq = lam * (pw * tp + mw * tn)
    local d1 = pw * e1p + mw * e1n
    local d2 = pw * e2p + mw * e2n
    local dd = d1 - d2
    sxx = tq + 2 * mu * (d2 + dd * 0.5 * (1 + c2))
    syy = tq + 2 * mu * (d2 + dd * 0.5 * (1 - c2))
    sxy = mu * dd * s2
    -- heaviside steps, 1 for positive values, 0 otherwise
    local ht = min(max(tr * 1e300, 0), 1)
    local h1 = min(max(e1 * 1e300, 0), 1)
    local h2 = min(max(e2 * 1e300, 0), 1)
    -- divided difference (e1+ - e2+)/(e1 - e2), tends to h1 for e1 = e2
    local eps = 1e-12 * (abs(e1) + abs(e2)) + 1e-300
    local qp = (e1p - e2p + h1 * eps) / (e1 - e2 + eps)
    local L = lam * (mw + (pw - mw) * ht)
    local K1 = 2 * mu * (mw + (pw - mw) * h1)
    local K2 = 2 * mu * (mw + (pw - mw) * h2)
    local K3 = mu * (pw * qp + mw * (1 - qp))
    local cc = 0.5 * (1 + c2)
    local ss = 0.5 * (1 - c2)
    local cs = 0.5 * s2
    local cs2 = cs * cs
    cg[o] = W * (L + K1 * cc * cc + K2 * ss * ss + 4 * K3 * cs2)
    cg[o + 1] = W * (L + K1 * ss * ss + K2 * cc * cc + 4 * K3 * cs2)
    cg[o + 2] = W * ((K1 + K2) * cs2 + K3 * c2 * c2)
    cg[o + 3] = W * (L + (K1 + K2) * cc * ss - 4 * K3 * cs2)
    cg[o + 4] = W * (K1 * cc * cs - K2 * ss * cs - 2 * K3 * cs * c2)
    cg[o + 5] = W * (K1 * ss * cs - K2 * cc * cs + 2 * K3 * cs * c2)
  else
    local tr = exx + eyy
    psip = 0.5 * lam * tr * tr + mu * (exx * exx + eyy * eyy + 0.5 * gxy * gxy)
    psim = 0
    sxx = pw * (lam * tr + 2 * mu * exx)
    syy = pw * (lam * tr + 2 * mu * eyy)
    sxy = pw * mu * gxy
    local c = W * pw
    cg[o] = c * (lam + 2 * mu)
    cg[o + 1] = c * (lam + 2 * mu)
    cg[o + 2] = c * mu
    cg[o + 3] = c * lam
    cg[o + 4] = 0
    cg[o + 5] = 0
  end
  return W * sxx, W * syy, W * sxy, psip, psim
end

-- energies only (line search, damage drive)
local function gp_psi(exx, eyy, gxy, lam, mu, split)
  if split then
    local tr = exx + eyy
    local hd = 0.5 * (exx - eyy)
    local exy = 0.5 * gxy
    local r = sqrt(hd * hd + exy * exy)
    local m = 0.5 * tr
    local e1 = m + r
    local e2 = m - r
    local tp = max(tr, 0)
    local tn = min(tr, 0)
    local e1p = max(e1, 0)
    local e2p = max(e2, 0)
    local e1n = min(e1, 0)
    local e2n = min(e2, 0)
    return 0.5 * lam * tp * tp + mu * (e1p * e1p + e2p * e2p),
      0.5 * lam * tn * tn + mu * (e1n * e1n + e2n * e2n)
  end
  local tr = exx + eyy
  return 0.5 * lam * tr * tr + mu * (exx * exx + eyy * eyy + 0.5 * gxy * gxy), 0
end

M.gp_split = gp_split
M.gp_psi = gp_psi

local function maybe_yield(S)
  local f = S.yield_fn
  if f then
    f(S)
  end
end

-- element kernels ----------------------------------------------------------
-- every kernel works on a list of elements (el, ne), a list of unknown nodes
-- (nl, nn) and the list of all nodes touched by those elements (rl, nr).
-- element e has the lower-left node e and nodes en2[e], en3[e], en4[e]
-- counter-clockwise. a global solve passes the full lists.
-- the loops run once per gauss point (q = 0..3) so that every loop body
-- stays small enough for the jit compiler.
-- gauss point q of an element sits at (s, t) in {g1, g2}^2:
--   q0 (g1, g1), q1 (g2, g1), q2 (g2, g2), q3 (g1, g2)
-- du/dx = cb*(u2-u1) + ct*(u3-u4), du/dy = cl*(u4-u1) + cr*(u3-u2)

local GP_CB = { G2, G2, G1, G1 }
local GP_CT = { G1, G1, G2, G2 }
local GP_CL = { G2, G1, G1, G2 }
local GP_CR = { G1, G2, G2, G1 }
local GP_N = {
  { Q22, Q12, Q11, Q12 },
  { Q12, Q22, Q12, Q11 },
  { Q11, Q12, Q22, Q12 },
  { Q12, Q11, Q12, Q22 },
}

local function assemble_pass(S, q, el, ne)
  local ux, uy, th, av = S.ux, S.uy, S.th, S.av
  local rx, ry, cg = S.rx, S.ry, S.cg
  local en2, en3, en4 = S.en2, S.en3, S.en4
  local lam, mu, split = S.lam, S.mu, S.split
  local ex0, ey0, eg0 = S.ex0, S.ey0, S.eg0
  local cb, ct, cl, cr = GP_CB[q + 1], GP_CT[q + 1], GP_CL[q + 1], GP_CR[q + 1]
  local w = GP_N[q + 1]
  local w1, w2, w3, w4 = w[1], w[2], w[3], w[4]
  local off = 6 * q * S.n
  local energy = 0
  for c = 0, ne - 1 do
    local n1 = el[c]
    local n2, n3, n4 = en2[n1], en3[n1], en4[n1]
    local x1, x2, x3, x4 = ux[n1], ux[n2], ux[n3], ux[n4]
    local y1, y2, y3, y4 = uy[n1], uy[n2], uy[n3], uy[n4]
    local dxx = cb * (x2 - x1) + ct * (x3 - x4)
    local dyx = cl * (x4 - x1) + cr * (x3 - x2)
    local dxy = cb * (y2 - y1) + ct * (y3 - y4)
    local dyy = cl * (y4 - y1) + cr * (y3 - y2)
    local pw = w1 * av[n1] + w2 * av[n2] + w3 * av[n3] + w4 * av[n4]
    local mw = w1 * th[n1] + w2 * th[n2] + w3 * th[n3] + w4 * th[n4]
    local sx, sy, sg, pp, pm = gp_split(cg, off + 6 * n1,
      dxx + ex0, dyy + ey0, dyx + dxy + eg0, pw, mw, lam, mu, split)
    energy = energy + W * (pw * pp + mw * pm)
    local fbx, ftx, flx, fqx = cb * sx, ct * sx, cl * sg, cr * sg
    local fby, fty, fly, fqy = cb * sg, ct * sg, cl * sy, cr * sy
    rx[n1] = rx[n1] - fbx - flx
    rx[n2] = rx[n2] + fbx - fqx
    rx[n3] = rx[n3] + ftx + fqx
    rx[n4] = rx[n4] - ftx + flx
    ry[n1] = ry[n1] - fby - fly
    ry[n2] = ry[n2] + fby - fqy
    ry[n3] = ry[n3] + fty + fqy
    ry[n4] = ry[n4] - fty + fly
  end
  return energy
end

-- energy gradient into rx, ry (valid on the unknown nodes), tangent into cg;
-- returns the energy of the listed elements plus the foundation energy of
-- the unknown nodes (the only part that a change of those nodes can alter)
local function assemble(S, el, ne, nl, nn, rl, nr)
  local ux, uy, al, th, av = S.ux, S.uy, S.al, S.th, S.av
  local rx, ry, kf = S.rx, S.ry, S.kfo
  local eta = S.eta
  for c = 0, nr - 1 do
    local i = rl[c]
    local b = 1 - al[i]
    av[i] = th[i] * ((1 - eta) * b * b + eta)
    rx[i] = 0
    ry[i] = 0
  end
  local energy = 0
  for c = 0, nn - 1 do
    local i = nl[c]
    local x, y, k = ux[i], uy[i], kf[i]
    rx[i] = k * x
    ry[i] = k * y
    energy = energy + 0.5 * k * (x * x + y * y)
  end
  for q = 0, 3 do
    energy = energy + assemble_pass(S, q, el, ne)
  end
  return energy
end

-- inverse diagonal of the tangent for jacobi preconditioning
-- derivative tables: for node k (0..3) and gauss index (0 -> g1, 1 -> g2)
local DXT = ffi.new("double[8]", { -G2, -G1, G2, G1, G1, G2, -G1, -G2 })
local DYT = ffi.new("double[8]", { -G2, -G1, -G1, -G2, G1, G2, G2, G1 })
-- gauss point q -> (t index, s index)
local GT = ffi.new("int[4]", { 0, 0, 1, 1 })
local GS = ffi.new("int[4]", { 0, 1, 1, 0 })

local function element_node(S, e, k)
  if k == 0 then
    return e
  elseif k == 1 then
    return S.en2[e]
  elseif k == 2 then
    return S.en3[e]
  end
  return S.en4[e]
end

local function compute_diag(S, el, ne, nl, nn, rl, nr)
  local cg, kf, mx, my = S.cg, S.kfo, S.mx, S.my
  for c = 0, nr - 1 do
    local i = rl[c]
    mx[i] = 0
    my[i] = 0
  end
  for c = 0, nn - 1 do
    local i = nl[c]
    mx[i] = kf[i]
    my[i] = kf[i]
  end
  local n6 = 6 * S.n
  for c = 0, ne - 1 do
    local e = el[c]
    local o = 6 * e
    for k = 0, 3 do
      local dxx, dyy = 0, 0
      for q = 0, 3 do
        local b = o + q * n6
        local dX = DXT[2 * k + GT[q]]
        local dY = DYT[2 * k + GS[q]]
        local c13, c23, c33 = cg[b + 4], cg[b + 5], cg[b + 2]
        dxx = dxx + cg[b] * dX * dX + 2 * c13 * dX * dY + c33 * dY * dY
        dyy = dyy + cg[b + 1] * dY * dY + 2 * c23 * dX * dY + c33 * dX * dX
      end
      local nk = element_node(S, e, k)
      mx[nk] = mx[nk] + dxx
      my[nk] = my[nk] + dyy
    end
  end
  for c = 0, nn - 1 do
    local i = nl[c]
    mx[i] = 1 / mx[i]
    my[i] = 1 / my[i]
  end
end

-- all four gauss points in one pass: the matvec is the inner loop of the
-- solver and this body is small enough for the jit compiler
local function matvec_elements(S, xx, xy, yx, yy, el, ne)
  local cg = S.cg
  local en2, en3, en4 = S.en2, S.en3, S.en4
  local s1, s2, s3 = 6 * S.n, 12 * S.n, 18 * S.n
  for c = 0, ne - 1 do
    local n1 = el[c]
    local n2, n3, n4 = en2[n1], en3[n1], en4[n1]
    local x1, x2, x3, x4 = xx[n1], xx[n2], xx[n3], xx[n4]
    local y1, y2, y3, y4 = xy[n1], xy[n2], xy[n3], xy[n4]
    local bx, tx, lx, qx = x2 - x1, x3 - x4, x4 - x1, x3 - x2
    local by, ty, ly, qy = y2 - y1, y3 - y4, y4 - y1, y3 - y2
    local dx1x, dx2x = G2 * bx + G1 * tx, G1 * bx + G2 * tx
    local dy1x, dy2x = G2 * lx + G1 * qx, G1 * lx + G2 * qx
    local dx1y, dx2y = G2 * by + G1 * ty, G1 * by + G2 * ty
    local dy1y, dy2y = G2 * ly + G1 * qy, G1 * ly + G2 * qy
    local o = 6 * n1
    local e1, e2, e3
    e1, e2, e3 = dx1x, dy1y, dy1x + dx1y
    local s0x = cg[o] * e1 + cg[o + 3] * e2 + cg[o + 4] * e3
    local s0y = cg[o + 3] * e1 + cg[o + 1] * e2 + cg[o + 5] * e3
    local s0g = cg[o + 4] * e1 + cg[o + 5] * e2 + cg[o + 2] * e3
    local o1 = o + s1
    e1, e2, e3 = dx1x, dy2y, dy2x + dx1y
    local s1x = cg[o1] * e1 + cg[o1 + 3] * e2 + cg[o1 + 4] * e3
    local s1y = cg[o1 + 3] * e1 + cg[o1 + 1] * e2 + cg[o1 + 5] * e3
    local s1g = cg[o1 + 4] * e1 + cg[o1 + 5] * e2 + cg[o1 + 2] * e3
    local o2 = o + s2
    e1, e2, e3 = dx2x, dy2y, dy2x + dx2y
    local s2x = cg[o2] * e1 + cg[o2 + 3] * e2 + cg[o2 + 4] * e3
    local s2y = cg[o2 + 3] * e1 + cg[o2 + 1] * e2 + cg[o2 + 5] * e3
    local s2g = cg[o2 + 4] * e1 + cg[o2 + 5] * e2 + cg[o2 + 2] * e3
    local o3 = o + s3
    e1, e2, e3 = dx2x, dy1y, dy1x + dx2y
    local s3x = cg[o3] * e1 + cg[o3 + 3] * e2 + cg[o3 + 4] * e3
    local s3y = cg[o3 + 3] * e1 + cg[o3 + 1] * e2 + cg[o3 + 5] * e3
    local s3g = cg[o3 + 4] * e1 + cg[o3 + 5] * e2 + cg[o3 + 2] * e3
    local fdx1x, fdx2x = s0x + s1x, s2x + s3x
    local fdy1y, fdy2y = s0y + s3y, s1y + s2y
    local fdy1x, fdy2x = s0g + s3g, s1g + s2g
    local fdx1y, fdx2y = s0g + s1g, s2g + s3g
    local fbx, ftx = G2 * fdx1x + G1 * fdx2x, G1 * fdx1x + G2 * fdx2x
    local flx, fqx = G2 * fdy1x + G1 * fdy2x, G1 * fdy1x + G2 * fdy2x
    local fby, fty = G2 * fdx1y + G1 * fdx2y, G1 * fdx1y + G2 * fdx2y
    local fly, fqy = G2 * fdy1y + G1 * fdy2y, G1 * fdy1y + G2 * fdy2y
    yx[n1] = yx[n1] - fbx - flx
    yx[n2] = yx[n2] + fbx - fqx
    yx[n3] = yx[n3] + ftx + fqx
    yx[n4] = yx[n4] - ftx + flx
    yy[n1] = yy[n1] - fby - fly
    yy[n2] = yy[n2] + fby - fqy
    yy[n3] = yy[n3] + fty + fqy
    yy[n4] = yy[n4] - fty + fly
  end
end

-- tangent times vector, (yx, yy) = K (xx, xy); x must vanish on rl \ nl
local function matvec(S, xx, xy, yx, yy, el, ne, nl, nn, rl, nr)
  local kf = S.kfo
  for c = 0, nr - 1 do
    local i = rl[c]
    yx[i] = 0
    yy[i] = 0
  end
  for c = 0, nn - 1 do
    local i = nl[c]
    local k = kf[i]
    yx[i] = k * xx[i]
    yy[i] = k * xy[i]
  end
  matvec_elements(S, xx, xy, yx, yy, el, ne)
end

-- tension energy psi+ per gauss point into S.psi, returns the element energy
local function energy_pass(S, q, el, ne)
  local tx, ty = S.ux, S.uy
  local av, th, psi = S.av, S.th, S.psi
  local en2, en3, en4 = S.en2, S.en3, S.en4
  local lam, mu, split = S.lam, S.mu, S.split
  local ex0, ey0, eg0 = S.ex0, S.ey0, S.eg0
  local cb, ct, cl, cr = GP_CB[q + 1], GP_CT[q + 1], GP_CL[q + 1], GP_CR[q + 1]
  local w = GP_N[q + 1]
  local w1, w2, w3, w4 = w[1], w[2], w[3], w[4]
  local energy = 0
  for c = 0, ne - 1 do
    local n1 = el[c]
    local n2, n3, n4 = en2[n1], en3[n1], en4[n1]
    local x1, x2, x3, x4 = tx[n1], tx[n2], tx[n3], tx[n4]
    local y1, y2, y3, y4 = ty[n1], ty[n2], ty[n3], ty[n4]
    local dxx = cb * (x2 - x1) + ct * (x3 - x4)
    local dyx = cl * (x4 - x1) + cr * (x3 - x2)
    local dxy = cb * (y2 - y1) + ct * (y3 - y4)
    local dyy = cl * (y4 - y1) + cr * (y3 - y2)
    local pp, pm = gp_psi(dxx + ex0, dyy + ey0, dyx + dxy + eg0, lam, mu, split)
    local pw = w1 * av[n1] + w2 * av[n2] + w3 * av[n3] + w4 * av[n4]
    local mw = w1 * th[n1] + w2 * th[n2] + w3 * th[n3] + w4 * th[n4]
    energy = energy + W * (pw * pp + mw * pm)
    psi[4 * n1 + q] = pp
  end
  return energy
end

-- energy at the current u (same parts as assemble) and the nodal damage
-- drive drv[i] = sum over gauss points of W * N_i * psi+, rebuilt on rl
local function energy_drive(S, el, ne, nl, nn, rl, nr)
  local ux, uy, kf, drv, psi = S.ux, S.uy, S.kfo, S.drv, S.psi
  local ew, esw, es = S.ew, S.esw, S.es
  local energy = 0
  for c = 0, nn - 1 do
    local i = nl[c]
    local x, y = ux[i], uy[i]
    energy = energy + 0.5 * kf[i] * (x * x + y * y)
  end
  for q = 0, 3 do
    energy = energy + energy_pass(S, q, el, ne)
  end
  for c = 0, nr - 1 do
    local i = rl[c]
    local a = 4 * i
    local b = 4 * ew[i]
    local d = 4 * esw[i]
    local f = 4 * es[i]
    drv[i] = W * (Q22 * psi[a] + Q12 * psi[a + 1] + Q11 * psi[a + 2] + Q12 * psi[a + 3]
      + Q12 * psi[b] + Q22 * psi[b + 1] + Q12 * psi[b + 2] + Q11 * psi[b + 3]
      + Q11 * psi[d] + Q12 * psi[d + 1] + Q22 * psi[d + 2] + Q12 * psi[d + 3]
      + Q12 * psi[f] + Q11 * psi[f + 1] + Q12 * psi[f + 2] + Q22 * psi[f + 3])
  end
  return energy
end

-- jacobi preconditioned conjugate gradient on the unknown nodes:
-- tangent * du = -gradient, with du = 0 on rl \ nl
local function pcg(S, el, ne, nl, nn, rl, nr)
  local rx, ry, dx, dy = S.rx, S.ry, S.dux, S.duy
  local zx, zy, px, py, qx, qy = S.zx, S.zy, S.px, S.py, S.qx, S.qy
  local mx, my = S.mx, S.my
  for c = 0, nr - 1 do
    local i = rl[c]
    px[i] = 0
    py[i] = 0
  end
  local rr, rz = 0, 0
  for c = 0, nn - 1 do
    local i = nl[c]
    local a, b = -rx[i], -ry[i]
    rx[i] = a
    ry[i] = b
    dx[i] = 0
    dy[i] = 0
    local za, zb = a * mx[i], b * my[i]
    px[i] = za
    py[i] = zb
    rr = rr + a * a + b * b
    rz = rz + a * za + b * zb
  end
  local r0 = sqrt(rr)
  S.last_res = r0
  local tol = S.cg_rtol * r0
  if tol < S.cg_atol then
    tol = S.cg_atol
  end
  if r0 <= tol then
    return 0
  end
  local it = 0
  local done = false
  while it < S.cg_maxit and not done do
    it = it + 1
    matvec(S, px, py, qx, qy, el, ne, nl, nn, rl, nr)
    local pq = 0
    for c = 0, nn - 1 do
      local i = nl[c]
      pq = pq + px[i] * qx[i] + py[i] * qy[i]
    end
    local alpha = rz / pq
    local rzn = 0
    rr = 0
    for c = 0, nn - 1 do
      local i = nl[c]
      dx[i] = dx[i] + alpha * px[i]
      dy[i] = dy[i] + alpha * py[i]
      local a, b = rx[i] - alpha * qx[i], ry[i] - alpha * qy[i]
      rx[i] = a
      ry[i] = b
      rr = rr + a * a + b * b
      local za, zb = a * mx[i], b * my[i]
      zx[i] = za
      zy[i] = zb
      rzn = rzn + a * za + b * zb
    end
    if sqrt(rr) <= tol then
      done = true
    else
      local beta = rzn / rz
      rz = rzn
      for c = 0, nn - 1 do
        local i = nl[c]
        px[i] = zx[i] + beta * px[i]
        py[i] = zy[i] + beta * py[i]
      end
      if it % 4 == 0 then
        maybe_yield(S)
      end
    end
  end
  S.cg_total = S.cg_total + it
  return it
end

-- one newton step for u with backtracking line search; leaves S.drv
-- consistent with the accepted displacement
local function newton_u(S, el, ne, nl, nn, rl, nr)
  local e0 = assemble(S, el, ne, nl, nn, rl, nr)
  compute_diag(S, el, ne, nl, nn, rl, nr)
  pcg(S, el, ne, nl, nn, rl, nr)
  local ux, uy, dx, dy = S.ux, S.uy, S.dux, S.duy
  local tau = 1
  local accepted = false
  local k = 0
  S.last_energy = e0
  while k < 8 and not accepted do
    k = k + 1
    for c = 0, nn - 1 do
      local i = nl[c]
      ux[i] = ux[i] + tau * dx[i]
      uy[i] = uy[i] + tau * dy[i]
    end
    local e1 = energy_drive(S, el, ne, nl, nn, rl, nr)
    if e1 <= e0 + 1e-11 * abs(e0) then
      accepted = true
      S.last_energy = e1
    else
      for c = 0, nn - 1 do
        local i = nl[c]
        ux[i] = ux[i] - tau * dx[i]
        uy[i] = uy[i] - tau * dy[i]
      end
      tau = 0.5 * tau
    end
  end
  if not accepted then
    energy_drive(S, el, ne, nl, nn, rl, nr)
  end
  S.newton_total = S.newton_total + 1
end

-- damage problem -------------------------------------------------------------

-- nodes that may change: damaged nodes, nodes close to the damage threshold,
-- and a neighborhood around them wide enough for the whole AT1 profile
-- (seed mask dilated by a square of radius R with two running-window
-- passes, cost independent of R)
local function build_candidates(S)
  local N, n = S.N, S.n
  local al, drv, th, c1 = S.al, S.drv, S.th, S.c1
  local mark, mark2, cand = S.mark, S.mark2, S.cand
  local fac = 2 * (1 - S.eta)
  local R = min(S.cand_r, floor((N - 1) / 2))
  for i = 0, n - 1 do
    if al[i] > 0 or fac * th[i] * drv[i] > 0.9 * c1[i] then
      mark[i] = 1
    else
      mark[i] = 0
    end
  end
  -- along x: mark2 = 1 where a seed lies within R in the same row
  for j = 0, N - 1 do
    local row = j * N
    local cnt = 0
    for d = -R, R do
      cnt = cnt + mark[row + d % N]
    end
    for i = 0, N - 1 do
      mark2[row + i] = cnt > 0 and 1 or 0
      cnt = cnt - mark[row + (i - R) % N] + mark[row + (i + R + 1) % N]
    end
  end
  -- along y
  local nc = 0
  for i = 0, N - 1 do
    local cnt = 0
    for d = -R, R do
      cnt = cnt + mark2[(d % N) * N + i]
    end
    for j = 0, N - 1 do
      mark[j * N + i] = cnt > 0 and 1 or 0
      cnt = cnt - mark2[((j - R) % N) * N + i] + mark2[((j + R + 1) % N) * N + i]
    end
  end
  for i = 0, n - 1 do
    if mark[i] == 1 then
      cand[nc] = i
      nc = nc + 1
    end
    mark[i] = 0
  end
  return nc
end

local function pgs_sweep(S, cl, first, last, step)
  local N, n = S.N, S.n
  local al, lb, drv, th, c1, lap = S.al, S.lb, S.drv, S.th, S.c1, S.lap
  local fac = 2 * (1 - S.eta)
  for c = first, last, step do
    local i = cl[c]
    local x = i % N
    local yc = i - x
    local xm = x == 0 and N - 1 or x - 1
    local xp = x == N - 1 and 0 or x + 1
    local ym = yc == 0 and n - N or yc - N
    local yp = yc == n - N and 0 or yc + N
    local o = i * 9
    local s = lap[o + 1] * al[yc + xp] + lap[o + 2] * al[yc + xm]
      + lap[o + 3] * al[yp + x] + lap[o + 4] * al[ym + x]
      + lap[o + 5] * al[yp + xp] + lap[o + 6] * al[yp + xm]
      + lap[o + 7] * al[ym + xp] + lap[o + 8] * al[ym + xm]
    local d = fac * th[i] * drv[i]
    local v = (d - c1[i] - s) / (d + lap[o])
    local l = lb[i]
    if v < l then
      v = l
    elseif v > 1 then
      v = 1
    end
    al[i] = v
  end
end

-- projected gauss-seidel on the candidate list; returns the largest change
-- and marks the tiles where damage changed by more than the tolerance
local function damage_step(S, cl, nc)
  local al, ao = S.al, S.ao
  local N, T, ntx = S.N, S.tile, S.ntx
  local tchg = S.tchg
  for c = 0, nc - 1 do
    ao[c] = al[cl[c]]
  end
  for sweep = 1, S.pgs_sweeps do
    pgs_sweep(S, cl, 0, nc - 1, 1)
    pgs_sweep(S, cl, nc - 1, 0, -1)
  end
  -- over-relaxed alternate minimization: same fixed points, fewer iterations
  local omega = S.omega
  if omega ~= 1 then
    local lb = S.lb
    for c = 0, nc - 1 do
      local i = cl[c]
      local v = ao[c] + omega * (al[i] - ao[c])
      local l = lb[i]
      if v < l then
        v = l
      elseif v > 1 then
        v = 1
      end
      al[i] = v
    end
  end
  -- a point damaged beyond a_open is part of an open crack: a = 1 exactly,
  -- and it no longer holds on to the substrate. the core of the smeared band
  -- stands for a crack of zero width; if it stayed bonded it would pin the
  -- crack faces to the substrate and the band would widen as the opening grows
  local a_open = S.a_open
  if a_open < 1 then
    local lb, kf, kfo, eta = S.lb, S.kf, S.kfo, S.eta
    for c = 0, nc - 1 do
      local i = cl[c]
      if al[i] >= a_open then
        al[i] = 1
        lb[i] = 1
        kfo[i] = eta * kf[i]
      end
    end
  end
  local dmax = 0
  local tol = S.tol_alpha
  for c = 0, nc - 1 do
    local i = cl[c]
    local d = abs(al[i] - ao[c])
    if d > dmax then
      dmax = d
    end
    if d >= tol then
      local x = i % N
      local y = (i - x) / N
      tchg[floor(y / T) * ntx + floor(x / T)] = 1
    end
  end
  S.last_candidates = nc
  return dmax
end

-- active region for local solves: changed tiles dilated by the solve radius.
-- active elements: all elements of active tiles; unknown nodes: nodes whose
-- four elements are all active; touched nodes: all nodes of active elements.
-- returns the number of active tiles, or -1 when a global solve is better
local function build_active(S)
  local N, T, ntx, nty = S.N, S.tile, S.ntx, S.nty
  local nt = ntx * nty
  local tchg, tact = S.tchg, S.tact
  local rt = S.tile_radius
  for t = 0, nt - 1 do
    tact[t] = 0
  end
  for t = 0, nt - 1 do
    if tchg[t] == 1 then
      tchg[t] = 0
      local tx = t % ntx
      local ty = (t - tx) / ntx
      for dy = -rt, rt do
        local row = ((ty + dy) % nty) * ntx
        for dx = -rt, rt do
          tact[row + (tx + dx) % ntx] = 1
        end
      end
    end
  end
  local count = 0
  for t = 0, nt - 1 do
    count = count + tact[t]
  end
  if count == 0 then
    return 0
  end
  if count > S.local_frac * nt then
    return -1
  end
  local ael, anl, enl, mark = S.ael, S.anl, S.enl, S.mark
  local en2, en3, en4 = S.en2, S.en3, S.en4
  local nae, nan, nen = 0, 0, 0
  for t = 0, nt - 1 do
    if tact[t] == 1 then
      local tx = t % ntx
      local ty = (t - tx) / ntx
      local i0, j0 = tx * T, ty * T
      local i1, j1 = min(i0 + T, N) - 1, min(j0 + T, N) - 1
      for j = j0, j1 do
        local tym = floor(((j - 1) % N) / T) * ntx
        local tyc = ty * ntx
        for i = i0, i1 do
          local e = j * N + i
          ael[nae] = e
          nae = nae + 1
          local txm = floor(((i - 1) % N) / T)
          if tact[tyc + txm] == 1 and tact[tym + txm] == 1 and tact[tym + tx] == 1 then
            anl[nan] = e
            nan = nan + 1
          end
          local k1, k2, k3, k4 = e, en2[e], en3[e], en4[e]
          if mark[k1] == 0 then
            mark[k1] = 1
            enl[nen] = k1
            nen = nen + 1
          end
          if mark[k2] == 0 then
            mark[k2] = 1
            enl[nen] = k2
            nen = nen + 1
          end
          if mark[k3] == 0 then
            mark[k3] = 1
            enl[nen] = k3
            nen = nen + 1
          end
          if mark[k4] == 0 then
            mark[k4] = 1
            enl[nen] = k4
            nen = nen + 1
          end
        end
      end
    end
  end
  for c = 0, nen - 1 do
    mark[enl[c]] = 0
  end
  S.nae, S.nan, S.nen = nae, nan, nen
  return count
end

-- fracture energy of the current damage field
function M.fracture_energy(S)
  local N, n = S.N, S.n
  local al, c1, lap = S.al, S.c1, S.lap
  local e = 0
  for i = 0, n - 1 do
    local a = al[i]
    if a > 0 then
      local x = i % N
      local yc = i - x
      local xm = x == 0 and N - 1 or x - 1
      local xp = x == N - 1 and 0 or x + 1
      local ym = yc == 0 and n - N or yc - N
      local yp = yc == n - N and 0 or yc + N
      local o = i * 9
      local s = lap[o] * a + lap[o + 1] * al[yc + xp] + lap[o + 2] * al[yc + xm]
        + lap[o + 3] * al[yp + x] + lap[o + 4] * al[ym + x]
        + lap[o + 5] * al[yp + xp] + lap[o + 6] * al[yp + xm]
        + lap[o + 7] * al[ym + xp] + lap[o + 8] * al[ym + xm]
      e = e + c1[i] * a + 0.5 * a * s
    end
  end
  return e
end

-- setup -----------------------------------------------------------------------

local function build_laplacian(S)
  local N, n = S.N, S.n
  local gc, th, lap, c1 = S.gc_n, S.th, S.lap, S.c1
  local ell = S.ell
  -- bilinear finite element couplings: diagonal, edge and corner
  local wd, we, wc = 2 / 3, 1 / 6, 1 / 3
  local ce = new_array(n, "double")
  for j = 0, N - 1 do
    local jp = (j + 1) % N
    for i = 0, N - 1 do
      local ip = (i + 1) % N
      local a = gc[j * N + i] * th[j * N + i] + gc[j * N + ip] * th[j * N + ip]
        + gc[jp * N + ip] * th[jp * N + ip] + gc[jp * N + i] * th[jp * N + i]
      ce[j * N + i] = 0.75 * ell * 0.25 * a
    end
  end
  for j = 0, N - 1 do
    local jm = (j - 1) % N
    for i = 0, N - 1 do
      local im = (i - 1) % N
      local k = j * N + i
      local ne = ce[j * N + i]
      local nw = ce[j * N + im]
      local sw = ce[jm * N + im]
      local se = ce[jm * N + i]
      local o = k * 9
      lap[o] = wd * (ne + nw + sw + se)
      lap[o + 1] = -we * (ne + se)
      lap[o + 2] = -we * (nw + sw)
      lap[o + 3] = -we * (ne + nw)
      lap[o + 4] = -we * (se + sw)
      lap[o + 5] = -wc * ne
      lap[o + 6] = -wc * nw
      lap[o + 7] = -wc * se
      lap[o + 8] = -wc * sw
      c1[k] = 0.375 * gc[k] * th[k] / ell
    end
  end
end

-- pre-existing micro-cracks: fully broken short segments. all four nodes of
-- every cell the segment passes are broken, so the segment cuts a strip of
-- elements and is traction free like a real crack
local function place_flaws(S, r)
  local N = S.N
  local al, lb = S.al, S.lb
  for f = 1, S.flaws do
    local cx = rng_uniform(r) * N
    local cy = rng_uniform(r) * N
    local ang = rng_uniform(r) * pi
    local half = 0.5 * S.flaw_len
    local steps = ceil(S.flaw_len * 2) + 1
    for s = 0, steps do
      local t = -half + S.flaw_len * s / steps
      local x = floor(cx + t * cos(ang)) % N
      local y = floor(cy + t * sin(ang)) % N
      local xp, yp = (x + 1) % N, (y + 1) % N
      al[y * N + x], lb[y * N + x] = 1, 1
      al[y * N + xp], lb[y * N + xp] = 1, 1
      al[yp * N + x], lb[yp * N + x] = 1, 1
      al[yp * N + xp], lb[yp * N + xp] = 1, 1
    end
  end
end

-- open cracks: nodes damaged beyond a_open are fully broken and no longer
-- bonded to the substrate. rebuilt from the damage field, so damage that is
-- set from outside (flaws, tests, a saved state) is treated the same way
function M.sync_open(S)
  local al, lb, kf, kfo = S.al, S.lb, S.kf, S.kfo
  local a_open, eta = S.a_open, S.eta
  for i = 0, S.n - 1 do
    if al[i] >= a_open then
      al[i] = 1
      lb[i] = 1
      kfo[i] = eta * kf[i]
    else
      kfo[i] = kf[i]
    end
  end
end

function M.new(params)
  local S = {}
  for k, v in pairs(M.defaults) do
    S[k] = v
  end
  if params then
    for k, v in pairs(params) do
      S[k] = v
    end
  end
  local N = S.n
  S.N = N
  S.n = N * N
  local n = S.n
  local E = 1
  local nu = S.nu
  S.lam = E * nu / (1 - nu * nu)
  S.mu = E / (2 * (1 + nu))
  S.ebar = E / (1 - nu * nu)
  local ang = S.angle * pi / 180
  local c, s = cos(ang), sin(ang)
  S.sxx = S.s1 * c * c + S.s2 * s * s
  S.syy = S.s1 * s * s + S.s2 * c * c
  S.sxy = (S.s1 - S.s2) * c * s
  local names = { "ux", "uy", "tx", "al", "lb", "ao", "th", "gc_n", "kf", "av", "drv",
    "rx", "ry", "zx", "zy", "px", "py", "qx", "qy", "dux", "duy", "mx", "my", "c1", "kfo" }
  for _, name in ipairs(names) do
    S[name] = new_array(n, "double")
  end
  S.lap = new_array(9 * n, "double")
  S.cg = new_array(24 * n, "double")
  S.psi = new_array(4 * n, "double")
  S.mark = new_array(n, "uint8_t")
  S.mark2 = new_array(n, "uint8_t")
  local inames = { "cand", "all", "en2", "en3", "en4", "ew", "esw", "es", "ael", "anl", "enl" }
  for _, name in ipairs(inames) do
    S[name] = new_array(n, "int32_t")
  end
  -- connectivity (periodic)
  for j = 0, N - 1 do
    local jp, jm = (j + 1) % N, (j - 1) % N
    for i = 0, N - 1 do
      local ip, im = (i + 1) % N, (i - 1) % N
      local e = j * N + i
      S.all[e] = e
      S.en2[e] = j * N + ip
      S.en3[e] = jp * N + ip
      S.en4[e] = jp * N + i
      S.ew[e] = j * N + im
      S.esw[e] = jm * N + im
      S.es[e] = jm * N + i
    end
  end
  -- tiles for local solves
  S.ntx = ceil(N / S.tile)
  S.nty = S.ntx
  S.tchg = new_array(S.ntx * S.nty, "uint8_t")
  S.tact = new_array(S.ntx * S.nty, "uint8_t")
  local r = rng_new(S.seed)
  S.rng = r
  -- thickness factor (log-normal, mean 1)
  if S.h_std > 0 then
    gaussian_field(N, S.h_corr, r, S.th, S.tx)
    local mean = 0
    for i = 0, n - 1 do
      S.th[i] = exp(S.h_std * S.th[i])
      mean = mean + S.th[i]
    end
    mean = mean / n
    for i = 0, n - 1 do
      S.th[i] = S.th[i] / mean
    end
  else
    for i = 0, n - 1 do
      S.th[i] = 1
    end
  end
  if S.thickness_map then
    -- optional user map (array of n values, mean about 1)
    for i = 0, n - 1 do
      S.th[i] = S.thickness_map[i]
    end
  end
  -- fracture energy (log-normal with median gc)
  if S.gc_std > 0 then
    gaussian_field(N, S.gc_corr, r, S.gc_n, S.tx)
    for i = 0, n - 1 do
      S.gc_n[i] = S.gc * exp(S.gc_std * S.gc_n[i])
    end
  else
    for i = 0, n - 1 do
      S.gc_n[i] = S.gc
    end
  end
  local le_max = 0
  for i = 0, n - 1 do
    S.tx[i] = 0
    S.kf[i] = S.ebar / (S.le * S.le * S.th[i] ^ S.kexp)
    local l = sqrt(S.ebar * S.th[i] / S.kf[i])
    if l > le_max then
      le_max = l
    end
  end
  S.tile_radius = ceil(S.local_radius * le_max / S.tile)
  build_laplacian(S)
  S.cand_r = ceil(2 * S.ell) + 1
  place_flaws(S, r)
  M.sync_open(S)
  -- load unit: damage onset of the homogeneous film with median properties
  local pp = gp_psi(S.sxx, S.syy, 2 * S.sxy, S.lam, S.mu, S.split)
  if pp <= 0 then
    error("craquelure: shrink tensor gives no tension, nothing can crack")
  end
  S.pref = sqrt(3 * S.gc / (16 * S.ell * (1 - S.eta) * pp))
  -- unit of load paths: onset of the same film under equal shrinkage
  local piso = gp_psi(1, 1, 0, S.lam, S.mu, S.split)
  S.mref = sqrt(3 * S.gc / (16 * S.ell * (1 - S.eta) * piso))
  S.ex0, S.ey0, S.eg0 = 0, 0, 0
  S.m1, S.m2, S.mangle = 0, 0, 0
  S.progress = 0
  S.path_index = 0
  S.last_energy = 0
  -- scratch for the derived fields (render), never touched by the solver
  S.rw1 = new_array(n, "double")
  S.rw2 = new_array(n, "double")
  S.rw3 = new_array(n, "double")
  S.p = 0
  S.step = 0
  S.cg_total = 0
  S.newton_total = 0
  S.iter_total = 0
  S.local_total = 0
  S.last_res = 0
  S.last_candidates = 0
  S.nae, S.nan, S.nen = 0, 0, 0
  S.cg_atol = 1e-9 * S.pref * sqrt(n)
  S.history = {}
  S.done = false
  S.stop = false
  return S
end

-- quasi-static evolution -----------------------------------------------------

local function newton_global(S)
  local all, n = S.all, S.n
  newton_u(S, all, n, all, n, all, n)
end

local function newton_local(S)
  newton_u(S, S.ael, S.nae, S.anl, S.nan, S.enl, S.nen)
end

-- mismatch tensor (xx, yy, engineering shear xy) from principal values,
-- s1 along angle (degrees) and s2 across it
local function principal_to_tensor(s1, s2, angle)
  local a = angle * pi / 180
  local c, s = cos(a), sin(a)
  return s1 * c * c + s2 * s * s, s1 * s * s + s2 * c * c, 2 * (s1 - s2) * c * s
end

-- set the mismatch (eigenstrain) tensor. under proportional loading with
-- fixed damage the displacement scales with the load (the energy is
-- quadratic), so the old solution is scaled as a predictor
local function set_mismatch(S, exx, eyy, gxy)
  local ox, oy, og = S.ex0, S.ey0, S.eg0
  local oo = ox * ox + oy * oy + 0.5 * og * og
  local nn = exx * exx + eyy * eyy + 0.5 * gxy * gxy
  if oo > 0 then
    local f = (exx * ox + eyy * oy + 0.5 * gxy * og) / oo
    local rx, ry, rg = exx - f * ox, eyy - f * oy, gxy - f * og
    if f > 0 and rx * rx + ry * ry + 0.5 * rg * rg <= 1e-18 * nn then
      local ux, uy = S.ux, S.uy
      for i = 0, S.n - 1 do
        ux[i] = ux[i] * f
        uy[i] = uy[i] * f
      end
    end
  end
  S.ex0, S.ey0, S.eg0 = exx, eyy, gxy
  -- principal values in units of the isotropic onset, for display
  local m = 0.5 * (exx + eyy)
  local hd = 0.5 * (exx - eyy)
  local r = sqrt(hd * hd + 0.25 * gxy * gxy)
  S.m1 = (m + r) / S.mref
  S.m2 = (m - r) / S.mref
  S.mangle = 0.5 * math.atan2(gxy, exx - eyy) * 180 / pi
end

-- minimize the energy at the current mismatch under irreversibility
local function solve_step(S)
  local n = S.n
  local al, lb = S.al, S.lb
  for i = 0, n - 1 do
    lb[i] = al[i]
  end
  M.sync_open(S)
  local iters = 0
  local dmax = 1
  local cg0 = S.cg_total
  local mode_local = false
  local converged = false
  while not converged and iters < S.max_iter and not S.stop do
    iters = iters + 1
    if mode_local then
      newton_local(S)
      dmax = damage_step(S, S.enl, S.nen)
      S.local_total = S.local_total + 1
      if dmax < S.tol_alpha then
        mode_local = false
      elseif build_active(S) <= 0 then
        mode_local = false
      end
    else
      newton_global(S)
      local nc = build_candidates(S)
      dmax = damage_step(S, S.cand, nc)
      if dmax < S.tol_alpha then
        converged = true
      elseif build_active(S) > 0 then
        mode_local = true
      end
    end
    S.iter_in_step = iters
    maybe_yield(S)
  end
  S.step = S.step + 1
  S.iter_total = S.iter_total + iters
  local broken = 0
  for i = 0, n - 1 do
    if al[i] > 0.5 then
      broken = broken + 1
    end
  end
  local h = {
    load = S.p / S.pref,
    m1 = S.m1,
    m2 = S.m2,
    iters = iters,
    cg = S.cg_total - cg0,
    broken = broken / n,
    e_el = S.last_energy,
    e_fr = M.fracture_energy(S),
    dmax = dmax,
  }
  S.history[#S.history + 1] = h
  return h
end

-- one quasi-static step: mismatch = p * (shrink tensor of the parameters)
function M.load_step(S, p)
  set_mismatch(S, p * S.sxx, p * S.syy, 2 * p * S.sxy)
  S.p = p
  return solve_step(S)
end

-- one quasi-static step with a mismatch given by its principal values in
-- units of the isotropic onset mref (s1 along angle in degrees, s2 across)
function M.mismatch_step(S, s1, s2, angle)
  local exx, eyy, gxy = principal_to_tensor(s1 * S.mref, s2 * S.mref, angle or 0)
  set_mismatch(S, exx, eyy, gxy)
  return solve_step(S)
end

-- load path: S.path is a list of mismatch states { s1, s2, angle } in units
-- of the isotropic onset; s2 defaults to s1 and angle to 0. the mismatch
-- moves linearly from state to state in steps of at most load_step per
-- component, so any history of drying, aging and humidity can be followed
local function run_path(S)
  local path = S.path
  local st = {}
  for k = 1, #path do
    local e = path[k]
    local s1 = e[1]
    local s2 = e[2] or s1
    local xx, yy, xy = principal_to_tensor(s1, s2, e[3] or 0)
    st[k] = { xx, yy, xy }
  end
  local total = 0
  for k = 2, #st do
    local a, b = st[k - 1], st[k]
    total = total + max(abs(b[1] - a[1]), abs(b[2] - a[2]), abs(b[3] - a[3]))
  end
  if total <= 0 then
    total = 1
  end
  local done = 0
  local mref = S.mref
  S.progress = 0
  S.path_index = 1
  set_mismatch(S, st[1][1] * mref, st[1][2] * mref, st[1][3] * mref)
  solve_step(S)
  local k = 2
  while k <= #st and not S.stop do
    local a, b = st[k - 1], st[k]
    local len = max(abs(b[1] - a[1]), abs(b[2] - a[2]), abs(b[3] - a[3]))
    local steps = max(1, ceil(len / S.load_step - 1e-9))
    local j = 1
    S.path_index = k
    while j <= steps and not S.stop do
      local t = j / steps
      set_mismatch(S, (a[1] + t * (b[1] - a[1])) * mref, (a[2] + t * (b[2] - a[2])) * mref,
        (a[3] + t * (b[3] - a[3])) * mref)
      S.progress = (done + t * len) / total
      solve_step(S)
      j = j + 1
    end
    done = done + len
    k = k + 1
  end
end

function M.run(S)
  if S.path then
    run_path(S)
  else
    local k = 0
    local p = S.load_start
    local span = S.load_end - S.load_start
    while p <= S.load_end * (1 + 1e-9) and not S.stop do
      S.progress = span > 0 and (p - S.load_start) / span or 1
      M.load_step(S, p * S.pref)
      k = k + 1
      p = S.load_start + k * S.load_step
    end
  end
  S.progress = 1
  S.done = true
end

-- solve the elastic state at load p without damage evolution (tests)
function M.solve_elastic(S, p, newton_steps)
  set_mismatch(S, p * S.sxx, p * S.syy, 2 * p * S.sxy)
  S.p = p
  M.sync_open(S)
  for k = 1, newton_steps or 3 do
    newton_global(S)
  end
end

-- derived fields -----------------------------------------------------------

-- element-centre elastic strain and stress of element (i, j)
function M.element_state(S, i, j)
  local N = S.N
  local ip, jp = (i + 1) % N, (j + 1) % N
  local n1, n2, n3, n4 = j * N + i, j * N + ip, jp * N + ip, jp * N + i
  local ux, uy = S.ux, S.uy
  local exx = 0.5 * ((ux[n2] - ux[n1]) + (ux[n3] - ux[n4])) + S.ex0
  local eyy = 0.5 * ((uy[n4] - uy[n1]) + (uy[n3] - uy[n2])) + S.ey0
  local gxy = 0.5 * ((ux[n4] - ux[n1]) + (ux[n3] - ux[n2]))
    + 0.5 * ((uy[n2] - uy[n1]) + (uy[n3] - uy[n4])) + S.eg0
  local b1, b2, b3, b4 = 1 - S.al[n1], 1 - S.al[n2], 1 - S.al[n3], 1 - S.al[n4]
  local eta = S.eta
  local pw = 0.25 * (S.th[n1] * ((1 - eta) * b1 * b1 + eta) + S.th[n2] * ((1 - eta) * b2 * b2 + eta)
    + S.th[n3] * ((1 - eta) * b3 * b3 + eta) + S.th[n4] * ((1 - eta) * b4 * b4 + eta))
  local mw = 0.25 * (S.th[n1] + S.th[n2] + S.th[n3] + S.th[n4])
  local tmp = S.state_tmp
  if not tmp then
    tmp = ffi.new("double[6]")
    S.state_tmp = tmp
  end
  local sxx, syy, sxy = gp_split(tmp, 0, exx, eyy, gxy, pw, mw, S.lam, S.mu, S.split)
  return exx, eyy, gxy, sxx / W, syy / W, sxy / W
end

-- crack opening at every node near a crack, in cells, into S.rw1.
-- in a broken element almost all of the strain is opening of the crack, so
-- an element adds (1 - stiffness) * (largest principal strain minus that of
-- the uncracked film) times its width to the displacement jump across the
-- band. a running box as wide as the band sums these across the band.
-- returns the opening of a single straight crack in the uncracked film,
-- 2 * sigma1 * le / ebar (shear lag), as a reference
function M.opening(S)
  local N, n = S.N, S.n
  local ux, uy, al = S.ux, S.uy, S.al
  local en2, en3, en4 = S.en2, S.en3, S.en4
  local w, tmp, out = S.rw2, S.rw3, S.rw1
  local ex0, ey0, eg0 = S.ex0, S.ey0, S.eg0
  local a_open = S.a_open
  local hd0 = 0.5 * (ex0 - ey0)
  local e1ref = 0.5 * (ex0 + ey0) + sqrt(hd0 * hd0 + 0.25 * eg0 * eg0)
  for e = 0, n - 1 do
    local n2, n3, n4 = en2[e], en3[e], en4[e]
    local exx = 0.5 * ((ux[n2] - ux[e]) + (ux[n3] - ux[n4])) + ex0
    local eyy = 0.5 * ((uy[n4] - uy[e]) + (uy[n3] - uy[n2])) + ey0
    local gxy = 0.5 * ((ux[n4] - ux[e]) + (ux[n3] - ux[n2]))
      + 0.5 * ((uy[n2] - uy[e]) + (uy[n3] - uy[n4])) + eg0
    local hd = 0.5 * (exx - eyy)
    local e1 = 0.5 * (exx + eyy) + sqrt(hd * hd + 0.25 * gxy * gxy)
    local b1, b2, b3, b4 = 1 - al[e], 1 - al[n2], 1 - al[n3], 1 - al[n4]
    local g = 0.25 * (b1 * b1 + b2 * b2 + b3 * b3 + b4 * b4)
    -- only elements at an open crack count: diffuse damage is not a crack
    local bmin = min(b1, b2, b3, b4)
    if bmin <= 1 - a_open then
      w[e] = max(e1 - e1ref, 0) * (1 - g)
    else
      w[e] = 0
    end
  end
  -- elements with lower-left node in [x - R - 1, x + R] around node x
  local R = ceil(2 * S.ell)
  local side = 2 * R + 2
  for j = 0, N - 1 do
    local row = j * N
    local sum = 0
    for d = -R - 1, R do
      sum = sum + w[row + d % N]
    end
    for i = 0, N - 1 do
      tmp[row + i] = sum
      sum = sum - w[row + (i - R - 1) % N] + w[row + (i + R + 1) % N]
    end
  end
  for i = 0, N - 1 do
    local sum = 0
    for d = -R - 1, R do
      sum = sum + tmp[(d % N) * N + i]
    end
    for j = 0, N - 1 do
      out[j * N + i] = sum / side
      sum = sum - tmp[((j - R - 1) % N) * N + i] + tmp[((j + R + 1) % N) * N + i]
    end
  end
  local tmp6 = S.state_tmp
  if not tmp6 then
    tmp6 = ffi.new("double[6]")
    S.state_tmp = tmp6
  end
  local sxx, syy, sxy = gp_split(tmp6, 0, ex0, ey0, eg0, 1, 1, S.lam, S.mu, S.split)
  sxx, syy, sxy = sxx / W, syy / W, sxy / W
  local s1 = 0.5 * (sxx + syy) + sqrt(0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy)
  return 2 * max(s1, 0) * S.le / S.ebar
end

-- islands: connected regions of intact film (damage below 0.5, periodic,
-- 4 neighbours). returns the number of islands with at least min_area nodes
-- and their mean area in cells
local function island_push(al, lab, stack, top, k, id)
  if lab[k] == 0 and al[k] < 0.5 then
    lab[k] = id
    stack[top] = k
    return top + 1
  end
  return top
end

function M.islands(S, min_area)
  local N, n = S.N, S.n
  local al = S.al
  min_area = min_area or 8
  local lab, stack = S.isl_lab, S.isl_stack
  if not lab then
    lab = new_array(n, "int32_t")
    stack = new_array(n, "int32_t")
    S.isl_lab, S.isl_stack = lab, stack
  end
  for i = 0, n - 1 do
    lab[i] = 0
  end
  local id, count, total = 0, 0, 0
  for s = 0, n - 1 do
    if lab[s] == 0 and al[s] < 0.5 then
      id = id + 1
      lab[s] = id
      stack[0] = s
      local top = 1
      local area = 0
      while top > 0 do
        top = top - 1
        local i = stack[top]
        area = area + 1
        local x = i % N
        local yc = i - x
        local xm = x == 0 and N - 1 or x - 1
        local xp = x == N - 1 and 0 or x + 1
        local ym = yc == 0 and n - N or yc - N
        local yp = yc == n - N and 0 or yc + N
        top = island_push(al, lab, stack, top, yc + xm, id)
        top = island_push(al, lab, stack, top, yc + xp, id)
        top = island_push(al, lab, stack, top, ym + x, id)
        top = island_push(al, lab, stack, top, yp + x, id)
      end
      if area >= min_area then
        count = count + 1
        total = total + area
      end
    end
  end
  return count, count > 0 and total / count or 0
end

-- render into an rgba8 buffer of (n*scale)^2 * 4 bytes, bilinear between
-- nodes (periodic). mode: "craquelure", "damage", "stress", "thickness".
-- in "craquelure" the width of a crack line grows with its opening,
-- width_scale sets how strongly (default 0.5, 0 = all lines equal)
function M.render(S, ptr, mode, scale, width_scale)
  local N, n = S.N, S.n
  scale = scale or 1
  width_scale = width_scale or 0.5
  local Wd = N * scale
  local al, th = S.al, S.th
  local buf = ffi.cast("uint8_t*", ptr)
  local f = S.tx
  local op = S.rw1
  local oinv = 0
  -- node field to show
  if mode == "stress" then
    -- largest principal stress of the element to the upper right of a node,
    -- in units of the largest principal stress of the uncracked film
    local tmp6 = S.state_tmp
    if not tmp6 then
      tmp6 = ffi.new("double[6]")
      S.state_tmp = tmp6
    end
    local sxx, syy, sxy = gp_split(tmp6, 0, S.ex0, S.ey0, S.eg0, 1, 1, S.lam, S.mu, S.split)
    sxx, syy, sxy = sxx / W, syy / W, sxy / W
    local ref = 0.5 * (sxx + syy) + sqrt(0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy)
    if ref <= 1e-300 then
      ref = 1
    end
    for j = 0, N - 1 do
      for i = 0, N - 1 do
        local _, _, _, ex, ey, ez = M.element_state(S, i, j)
        local m = 0.5 * (ex + ey)
        local r = sqrt(0.25 * (ex - ey) * (ex - ey) + ez * ez)
        f[j * N + i] = (m + r) / ref
      end
    end
  elseif mode == "thickness" then
    local lo, hi = 1e300, -1e300
    for i = 0, n - 1 do
      lo = min(lo, th[i])
      hi = max(hi, th[i])
    end
    local span = hi - lo
    if span < 1e-12 then
      span = 1
    end
    for i = 0, n - 1 do
      f[i] = (th[i] - lo) / span
    end
  else
    for i = 0, n - 1 do
      f[i] = al[i]
    end
    if mode == "craquelure" then
      local oref = M.opening(S)
      if oref > 1e-300 then
        oinv = 1 / oref
      end
    end
  end
  local inv = 1 / scale
  for py = 0, Wd - 1 do
    local y = (py + 0.5) * inv - 0.5
    local j0 = floor(y)
    local fy = y - j0
    local r0 = (j0 % N) * N
    local r1 = ((j0 + 1) % N) * N
    for px = 0, Wd - 1 do
      local x = (px + 0.5) * inv - 0.5
      local i0 = floor(x)
      local fx = x - i0
      local c0 = i0 % N
      local c1 = (i0 + 1) % N
      local w00, w10 = (1 - fy) * (1 - fx), (1 - fy) * fx
      local w01, w11 = fy * (1 - fx), fy * fx
      local v = w00 * f[r0 + c0] + w10 * f[r0 + c1] + w01 * f[r1 + c0] + w11 * f[r1 + c1]
      local rr, gg, bb
      if mode == "damage" then
        local g = 255 * (1 - v)
        rr, gg, bb = g, g, g
      elseif mode == "thickness" then
        local g = 40 + 200 * v
        rr, gg, bb = g, g, g
      elseif mode == "stress" then
        if v >= 0 then
          v = min(v, 1)
          rr, gg, bb = 255, 255 * (1 - v), 255 * (1 - v)
        else
          v = min(-v, 1)
          rr, gg, bb = 255 * (1 - v), 255 * (1 - v), 255
        end
      else
        -- craquelure: dark line along the damage band, wider where the
        -- crack is open wider (old cracks in large islands), a hairline
        -- where it is closed or barely open. the paint tone follows the
        -- local film thickness
        local q = (w00 * op[r0 + c0] + w10 * op[r0 + c1] + w01 * op[r1 + c0] + w11 * op[r1 + c1]) * oinv
        q = min(max(q, 0), 1.5)
        local vt = min(max(0.95 - 0.8 * width_scale * q, 0.1), 0.95)
        local t = (v - vt) / 0.2 + 0.5
        t = min(max(t, 0), 1)
        local d = t * t * (3 - 2 * t)
        local h = w00 * th[r0 + c0] + w10 * th[r0 + c1] + w01 * th[r1 + c0] + w11 * th[r1 + c1]
        local tone = 0.94 + 0.06 * (h - 1)
        rr = (214 * tone) * (1 - d) + 46 * d
        gg = (182 * tone) * (1 - d) + 30 * d
        bb = (122 * tone) * (1 - d) + 18 * d
      end
      local o = 4 * (py * Wd + px)
      buf[o] = min(max(rr, 0), 255)
      buf[o + 1] = min(max(gg, 0), 255)
      buf[o + 2] = min(max(bb, 0), 255)
      buf[o + 3] = 255
    end
  end
end

-- exposed for tests
M.set_mismatch = set_mismatch
M.assemble = assemble
M.energy_drive = energy_drive
M.compute_diag = compute_diag
M.matvec = matvec
M.damage_step = damage_step
M.build_candidates = build_candidates
M.newton_u = newton_u
M.newton_global = newton_global
M.gaussian_field = gaussian_field
M.build_laplacian = build_laplacian
M.new_array = new_array

return M
