function love.conf(t)
  t.identity = "craquelure"
  t.version = "11.4"
  t.window.title = "craquelure - phase-field film fracture"
  t.window.width = 1100
  t.window.height = 800
  t.window.resizable = true
  t.window.vsync = 1
  t.modules.joystick = false
  t.modules.physics = false
  t.modules.audio = false
  t.modules.sound = false
  t.modules.video = false
end
