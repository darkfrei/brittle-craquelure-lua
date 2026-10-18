-- presets for the craquelure model, all lengths in grid cells
-- le   shear-lag length of the unit film: crack spacing scales with it,
--      le = Z * thickness, Z ~ 2 on a substrate as stiff as the film,
--      larger on a softer one
-- ell  AT1 material length, sets the strength sigma_c = sqrt(3*gc*E/(8*ell))
-- le/ell  brittleness of the film relative to its thickness:
--      about 4: a crack can only grow when the load rises, many cracks
--               start together and grow toward each other
--      about 8: a crack that starts runs through the whole film at once,
--               the pattern forms crack after crack
-- s1, s2, angle  principal shrinkage mismatch (film minus substrate),
--      negative = the film is compressed in that direction; with
--      load_start, load_end in units of p_ref, the mismatch at which a
--      flawless homogeneous film under this mismatch starts to damage
-- path  instead of s1, s2 and a load range: a list of mismatch states
--      { along x, along y } in units of the isotropic onset (equal
--      shrinkage at which a flawless film starts to damage); the
--      mismatch moves from state to state in steps of load_step
-- wood panels: the grain runs along x. wood hardly changes its length
-- along the grain but swells and shrinks strongly across it

return {
  {
    name = "drying paint layer",
    info = {
      "equal shrinkage in all directions, le/ell = 4",
      "the film reaches its strength before a crack",
      "can run: many cracks start together, grow",
      "toward each other and meet at wide angles",
    },
    params = {
      ell = 2.0, le = 8.0, nu = 0.3,
      s1 = 1.0, s2 = 1.0, angle = 0,
      gc_std = 0.15, gc_corr = 3.0,
      load_start = 0.9, load_end = 1.5, load_step = 0.01,
    },
  },
  {
    name = "brittle glaze crazing",
    info = {
      "equal shrinkage, le/ell = 8 (thick brittle layer)",
      "a crack that starts runs across the film;",
      "later cracks split the islands and stop",
      "at older cracks (T junctions)",
    },
    params = {
      ell = 1.5, le = 12.0, nu = 0.25,
      s1 = 1.0, s2 = 1.0, angle = 0,
      gc_std = 0.1, gc_corr = 3.0,
      load_start = 0.85, load_end = 1.4, load_step = 0.01,
    },
  },
  {
    name = "panel painting, dry air",
    info = {
      "the gesso shrinks; along the grain (x) the wood",
      "keeps its length and stretches it, across the",
      "grain the wood shrinks more and compresses it:",
      "cracks run across the grain (Italian panels)",
    },
    params = {
      ell = 1.5, le = 12.0, nu = 0.3,
      s1 = 1.0, s2 = -0.5, angle = 0,
      gc_std = 0.1, gc_corr = 3.0,
      load_start = 0.85, load_end = 1.4, load_step = 0.01,
    },
  },
  {
    name = "panel painting, humid air",
    info = {
      "the dried gesso is stretched a little in all",
      "directions; then the wood swells across the",
      "grain (y) and stretches it more: cracks run",
      "along the grain (Flemish panels)",
    },
    params = {
      ell = 1.5, le = 12.0, nu = 0.3,
      gc_std = 0.1, gc_corr = 3.0,
      path = {
        { 0.85, 0.85 },
        { 0.75, 1.90 },
      },
      load_step = 0.01,
    },
  },
  {
    name = "panel painting, seasons",
    info = {
      "humid and dry seasons alternate: humid air",
      "opens cracks along the grain, dry air across",
      "it; closed cracks still carry compression.",
      "together they cut the gesso into a grid",
    },
    params = {
      ell = 1.5, le = 12.0, nu = 0.3,
      gc_std = 0.1, gc_corr = 3.0,
      path = {
        { 0.80, 0.80 },
        { 0.70, 1.80 },
        { 0.80, 0.80 },
        { 1.70, 0.10 },
        { 0.80, 0.80 },
        { 0.60, 2.20 },
        { 0.80, 0.80 },
        { 2.00, -0.10 },
      },
      load_step = 0.02,
    },
  },
  {
    name = "canvas painting, uneven tension",
    info = {
      "warp and weft pull differently:",
      "mismatch 1.0 along x and 0.7 along y",
      "first cracks cross the x direction,",
      "later ones connect them (Dutch canvases)",
    },
    params = {
      ell = 1.5, le = 10.0, nu = 0.3,
      s1 = 1.0, s2 = 0.7, angle = 0,
      gc_std = 0.1, gc_corr = 3.0,
      load_start = 0.85, load_end = 1.5, load_step = 0.01,
    },
  },
  {
    name = "impasto, thickness varies",
    info = {
      "film thickness varies about 2.5 fold",
      "le grows with thickness: thick areas break",
      "into large islands with wide cracks, thin",
      "areas into small ones ('v': thickness map)",
    },
    params = {
      ell = 1.5, le = 12.0, nu = 0.3,
      s1 = 1.0, s2 = 1.0, angle = 0,
      gc_std = 0.1, gc_corr = 3.0,
      h_std = 0.25, h_corr = 40.0, kexp = 1.0,
      load_start = 0.8, load_end = 1.5, load_step = 0.01,
    },
  },
  {
    name = "film with micro-cracks",
    info = {
      "60 short pre-existing cracks (flaws)",
      "cracks start from flaws below the strength of",
      "the intact film, as soon as the energy release",
      "of a crack channelling from a flaw reaches gc",
    },
    params = {
      ell = 2.0, le = 10.0, nu = 0.3,
      s1 = 1.0, s2 = 1.0, angle = 0,
      gc_std = 0.05, gc_corr = 3.0,
      flaws = 60, flaw_len = 10.0,
      load_start = 0.6, load_end = 1.4, load_step = 0.01,
    },
  },
}
