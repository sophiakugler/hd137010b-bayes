#import "@preview/physica:0.9.8": *
#import "@preview/frame-it:1.1.2": *
#import "@preview/big-todo:0.2.0": *
#import "@preview/unify:0.7.1": num, qty


#let e0 = $epsilon_0$
#let a0 = $a_0$
#let m0 = $mu_0$
#let anabla = $arrow(nabla)$
#let ea = $epsilon.alt$
#let ep = $epsilon$
#let m1 = $-1$
#let id = $bb(1)$
#let circ = $compose$
#let mattyp(x, y) = $"Mat"_(#x times #y)$
#let aut(x) = $"Aut"_(#x)$
#let adj(x) = $#x^\#$
#let prod = $times$
#let prodr = $times.r$
#let prodl = $times.l$
#let oprod = $times.o$
#let oadd = $plus.o$
#let iso = $tilde.equiv$
#let tensor = $oprod$


#let upar = $arrow.t$
#let doar = $arrow.b$

#let hspace = $#h(10pt)$
#let hsmall = $#h(5pt)$
#let hfill = box(width: 1fr, move(dy: -0.25em, line(length: 100%, stroke: 0.13pt + black)))

#let aJ = $arrow(J)$
#let aE = $arrow(E)$
#let aB = $arrow(B)$
#let ax = $arrow(x)$
#let ay = $arrow(y)$
#let az = $arrow(z)$
#let aF = $arrow(F)$
#let aa = $arrow(a)$
#let ab = $arrow(b)$
#let ac = $arrow(c)$
#let av = $arrow(v)$
#let ar = $arrow(r)$
#let aS = $arrow(S)$
#let ad = $arrow(d)$
#let aR = $arrow(R)$
#let ap = $arrow(p)$
#let an = $arrow(n)$
#let ak = $arrow(k)$
#let al = $arrow(l)$
#let au = $arrow(u)$
#let aA = $arrow(A)$
#let ass = $arrow(s)$
#let amu = $arrow(mu)$
#let am = $arrow(m)$
#let aM = $arrow(M)$
#let af = $arrow(f)$
#let adv = $dot(arrow(v))$
#let addx = $dot.double(arrow(x))$
#let adx = $dot(arrow(x))$
#let aPhi = $arrow(Phi)$
#let aP = $arrow(P)$
#let aH = $arrow(H)$
#let aD = $arrow(D)$
#let aB = $arrow(B)$
#let aC = $arrow(C)$
#let aphi = $arrow(phi)$
#let atheta = $arrow(theta)$
#let aK = $arrow(K)$
#let aj = $arrow(j)$
#let ai = $arrow(i)$
#let am = $arrow(m)$
#let aG = $arrow(G)$
#let ag = $arrow(g)$
#let aq = $arrow(q)$
#let aT = $arrow(T)$
#let aL = $arrow(L)$


#let hH = $hat(H)$
#let hA = $hat(A)$
#let hB = $hat(B)$
#let hP = $hat(P)$
#let hX = $hat(X)$
#let hE = $hat(E)$
#let hT = $hat(T)$
#let Psic = $Psi^*$


#let SO = $S O$
#let SU = $S U$
#let SO = $S O$
#let Sp = $S p$
#let SG = $S G$
#let SL = $S L$
#let GL = $G L$
#let IO = $I O$
#let io = $i o$
#let gl = $g l$
#let su = $s u$
#let so = $s o$
#let sp = $s p$
#let Ad = $A d$
#let ad = $a d$
#let sl = $s l$

#let Hom = $"Hom"$
#let Ker = $"Ker"$
#let Im = $"Im"$
#let Stab = $"Stab"$
#let tr = $"tr"$
#let span = $"span"$
#let Bij = $"Bij"$
#let End = $"End"$
#let rk = $"rk"$
#let Lie = $"Lie"$
#let Orb = $"Orb"$
#let Fix = $"Fix"$
#let Mat = $"Mat"$

#let ein = $wj$


#let arr(x) = $arrow(#x)$
#let adt(x) = $dot(arrow(#x))$
#let add(x) = $dot.double(arrow(#x))$
#let ddot(x) = $dot.double(#x)$
#let eval(x) = $evaluated(#x)$

#let ddot = $dot dot$
#let dddot = $dot dot dot$

#let avg(g) = $chevron.l #g chevron.r$

// Thermodynamic partial derivative $(diff^n x \/ diff y^n)_z$, held at constant
// `z`. `degree` is the optional 4th positional argument (an integer, default
// 1) -- same meaning as physica's own `pdv(f, x, degree)`.
//   #tpdv(S, T, V)     -> (diff S \/ diff T)_V
//   #tpdv(S, T, V, 2)  -> (diff^2 S \/ diff T^2)_V
#let tpdv(x, y, z, ..rest) = {
  let degree = rest.pos().at(0, default: 1)
  $(pdv(#x, #y, #degree))_#z$
}

#let ds2 = $dd(s, 2)$
#let ds = $dd(s)$
#let dt = $dd(t)$
#let dx = $dd(x)$
#let dx3 = $dd(x, 3)$
#let dy = $dd(y)$
#let dz = $dd(z)$

#let invhat(f) = $caron(#f)$
#let over(f) = $overline(#f)$
#let under(f) = $underline(#f)$
#let wave = $square$
#let waved = $square_delta$
#let cross = $crossproduct$
#let deg = $degree$

#let kg = $"kg"$
#let pc = $"pc"$
#let ly = $"ly"$
#let Mpc = $"Mpc"$
#let sr = $"sr"$
#let Hz = $"Hz"$
#let eV = $"eV"$
#let keV = $"keV"$
#let MeV = $"MeV"$
#let GeV = $"GeV"$
#let TeV = $"TeV"$
#let mm = $"mm"$
#let cm = $"cm"$
#let km = $"km"$
#let nm = $"nm"$
#let pm = $"pm"$

#let gloss = ()


// ---------------------------------------------------------------------------
// Breakouts (escape the column layout)
// ---------------------------------------------------------------------------
//
// `#colspan(n)[...]` emits a single invisible marker carrying its body.
// `cheatsheet-setup` reacts to it differently depending on `break_sections`:
//
//   break_sections: false (default) -- the custom packer spots the marker,
//   cuts the document there and lays the block out over `n` columns' worth of
//   width (n columns plus the gutters between them): the columns stop right
//   before it and restart right after it. Which columns the strip lands on is
//   worked out during layout: the boxes before it are balanced over the
//   columns it will cover, and the strip goes under them, as far left as it
//   can while still covering the column the text ran out in. A `#colspan(2)`
//   of three columns therefore sits over columns 2-3 when the text ran out in
//   column 3. The columns it does not cover keep filling beside it.
//
//   break_sections: true -- there is no packer pass (Typst's own column flow
//   lays the document out), so the marker is instead turned into a float
//   (`scope: "parent"`, same mechanism as `fullwidth-float` below) that spans
//   `n` columns' worth of width, always measured from the left. This is
//   content-blind: it always sinks to the top of whatever page it lands on,
//   in source order, rather than sitting where the text ran out. Reading the
//   marker's own resolved position (`it.location()`) to do better was tried
//   and reverted -- see the long comment in `cheatsheet-setup` where this
//   marker is handled. Short version: it is not safe. On a real multi-section
//   document, a debug dump of the computed placement and the placement Typst
//   actually rendered disagreed, meaning `it.location()` inside a show rule
//   that itself emits a float does not reliably converge to a value
//   consistent with the final layout. A float also always reserves its full
//   height across *every* column for that vertical band, even the ones it
//   does not cover. Getting the packer's precision needs the packer's full
//   control over layout, which `break_sections: true` trades away.
//
//   #colspan(2)[#image("images/wide-table.png", width: 100%)]
//   #fullwidth[...]                 // shorthand for every column
//
// `n` is clamped to the document's column count, so `#colspan(99)[...]` is
// full width. A breakout has to sit at the top level of the document -- one
// buried inside e.g. `#section[...]` is not seen by the splitter and simply
// stays in its column (in either mode).
//
// Named args:
//   gap    -- vertical space above/below the breakout (length or `auto`)

#let _fullwidth_key = "cheatsheet:fullwidth"

#let colspan(cols, body, gap: auto) = [
  #metadata((kind: _fullwidth_key, cols: cols, gap: gap, align: left, body: body))
  <cheatsheet-breakout>
]

// Every column, i.e. the whole text width. Nothing can sit beside it, so here
// `align` still means something: where a body narrower than the page sits.
#let fullwidth(body, gap: auto, align: center) = [
  #metadata((kind: _fullwidth_key, cols: auto, gap: gap, align: align, body: body))
  <cheatsheet-breakout>
]

#let _is_fullwidth(c) = (
  c.func() == metadata
    and type(c.value) == dictionary
    and c.value.at("kind", default: none) == _fullwidth_key
)

// Alternative that does *not* cut the columns: hand the body to Typst's float
// machinery so it spans all columns at the top / bottom of the page it lands
// on. Slower to place, but it can never collide with the column text.
//   #fullwidth-float(image("x.png"))            // top of page
//   #fullwidth-float(at: bottom, image("x.png"))
#let fullwidth-float(body, at: top, clearance: 0.6em) = place(
  at + center,
  float: true,
  scope: "parent",
  clearance: clearance,
  block(width: 100%, body),
)


// Horizontally tiles its arguments edge to edge, with no gap between them,
// scaling all of them (preserving each one's own aspect ratio) to a common
// height chosen so the row exactly fills `width` -- a justified photo-row
// layout. This is what keeps e.g. two #image calls from overflowing the
// column they sit in.
//   #hstack[#image("a.png")][#image("b.png")]
// `width` takes a fraction of the surrounding container's width.
#let hstack(..elements, spacing: 0pt, align: horizon, width: 100%) = layout(size => {
  let items = elements.pos()
  let dims = items.map(e => measure(e))
  let aspects = dims.map(d => d.width / d.height)
  let avail = size.width * (width / 100%) - spacing * (items.len() - 1)
  let h = avail / aspects.sum()
  let scaled = items.zip(dims).map(((e, d)) => {
    let factor = h / d.height
    box(width: d.width * factor, height: h, scale(
      x: factor * 100%,
      y: factor * 100%,
      origin: top + left,
      reflow: true,
      e,
    ))
  })
  stack(dir: ltr, spacing: spacing, ..scaled.map(e => std.align(align, e)))
})

// Define colors for different section numbers (curated palette)
#let section_colors = (
  rgb("#FACFDB"), // 1  - rose
  rgb("#C3DDFB"), // 7  - cornflower
  rgb("#BCE3CE"), // 2  - jade
  rgb("#EDD1EF"), // 3  - orchid
  rgb("#F0D6B7"), // 17 - amber
  rgb("#D0DFBC"), // 4  - moss
  rgb("#D9D6FB"), // 5  - periwinkle
  rgb("#FCCFD0"), // 10 - blush
  rgb("#B4E2EE"), // 9  - sky
  rgb("#F7D3BD"), // 8  - apricot
  rgb("#B5E3D9"), // 11 - seafoam
  rgb("#F5D0E5"), // 12 - fuchsia
  rgb("#C5E1C4"), // 13 - spring green
  rgb("#E4D4F6"), // 14 - lilac
  rgb("#DCDCB7"), // 15 - olive
  rgb("#CDDAFD"), // 16 - hyacinth
  rgb("#BADFF6"), // 18 - azure
  rgb("#FBD1C5"), // 19 - shell
  rgb("#B2E3E4"), // 20 - turquoise
  rgb("#E7D9B5"), // 6  - straw
)

// WCAG relative luminance, and the contrast ratio white text would have on a
// given background. 4.5 is the readability threshold for body-sized text.
#let _luminance(col) = {
  let (r, g, b) = rgb(col).components(alpha: false)
  let f = v => {
    let x = v / 100%
    if x <= 0.04045 { x / 12.92 } else { calc.pow((x + 0.055) / 1.055, 2.4) }
  }
  0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)
}

#let _white_contrast(col) = 1.05 / (_luminance(col) + 0.05)

// How washed out a box fill is. `0%` leaves the palette entry alone; positive
// values mix it toward white (more pastel, `100%` is white); negative values
// push the chroma back up (deeper, more saturated). Two branches rather than one
// mix, because "less pastel" means more colour, not more black.
#let pastel_tint(base, amount) = {
  if amount == 0% { base } else if amount > 0% { rgb(base).lighten(amount) } else {
    rgb(base).saturate(-amount)
  }
}

// The palette above holds the *full strength* section hues. They are what the
// title bar and the border are derived from; the box fill is the entry
// lightened by this much, which is where the sheet's pastel look comes from.
// `section_pastel` then moves the fill relative to it, so 0% is the intended
// look and the knob still runs in both directions.
//
// Why the literals are not simply pre-lightened: `section_bar_color` reads its
// chroma off whatever it is handed, and a 50% lightened entry has roughly half
// the chroma, which drops it under the 0.13 cap. Every bar would come out
// washed out with it -- jade's #008859 turns into #508169, sky's #0081a0 into
// #437f8e. The bar has to key off the saturated entry to keep its punch.
#let _fill_pastel = 50%

// The title bar's colour: the same hue as the box, with the chroma pushed up so
// it still reads as that section's colour, and darkened exactly as far as it
// has to be for white text to clear `ratio`. Bisects on lightness rather than
// picking a fixed one, since how dark a hue has to get varies a lot -- yellow
// has to go much further down than blue.
#let section_bar_color(base, ratio: 4.5) = {
  let (l, chroma, hue, ..) = oklch(base).components()
  let c = calc.min(chroma * 2.6, 0.13)
  let lo = 15% // dark enough for white text at any hue
  let hi = l
  for _ in range(14) {
    let mid = (lo + hi) / 2
    if _white_contrast(oklch(mid, c, hue)) >= ratio { lo = mid } else { hi = mid }
  }
  oklch(lo, c, hue)
}

// The section colour depends on the current level 1 heading, which is only
// known during layout. Rather than wrapping the block in `context` -- which
// makes it an opaque lump that the column packer cannot measure -- the block is
// emitted plainly and labelled; `cheatsheet-setup` paints it via a
// `show <cheatsheet-section>` rule.
//
// A titled box gets its heading turned into a full-width bar in the darkened
// colour with white text, so the start of a section is easy to pick out. An
// untitled continuation box has no heading and so no bar.
// `gap` is the space kept above and below the box, `radius` its corner
// rounding; both default to zero so that boxes tile as plain rectangles.
//
// A stroke is centred on the edge it is drawn along, so with the boxes tiled
// edge to edge each border would spill half its width over the neighbour. The
// negative outset pulls the painted rectangle in by exactly that half, which
// leaves the ink flush with the box's own layout area: borders meet, none
// overlaps, and no sliver of page shows through.
//
// `border`, `background` and `title_color` switch the three coloured parts off
// independently, for a plain black on white print.
//
// `pastel` moves the box fill relative to `_fill_pastel`, positive toward white
// and negative toward the full strength hue. It never touches the bar or the
// border, so white bar text keeps its contrast however pale the box gets.
#let section_fill(
  body,
  title_size: 11.5pt,
  gap: 0pt,
  radius: 0pt,
  stroke_width: 0.4pt,
  border: true,
  background: true,
  title_color: true,
  pastel: 0%,
) = context {
  let current_section = counter(heading.where(level: 1)).get().first()
  let color = if current_section >= 1 and current_section <= section_colors.len() {
    section_colors.at(current_section - 1)
  } else {
    gray.lighten(80%) // fallback color
  }

  let bar = section_bar_color(color)
  show heading.where(level: 2): it => block(
    width: 100%,
    fill: if title_color { bar } else { none },
    radius: (top-left: radius, top-right: radius),
    inset: (x: 4pt, y: 2.5pt),
    above: 0pt,
    below: 0pt,
    text(fill: if title_color { white } else { black }, size: title_size, weight: "bold", {
      if it.numbering != none [#counter(heading).display(it.numbering) ]
      it.body
    }),
  )

  block(
    width: body.at("width", default: 100%),
    fill: if background { pastel_tint(color, _fill_pastel + pastel) } else { none },
    stroke: if border { stroke_width + bar } else { none },
    // Without a border there is no ink to pull back in, and shrinking the fill
    // would open a sliver between neighbouring boxes.
    outset: if border { -stroke_width / 2 } else { 0pt },
    radius: radius,
    above: gap,
    below: gap,
    breakable: true,
    body,
  )
}

// Main section function. Two forms, nothing else:
//
//   #section[Title][Content]  -- titled box, numbered off the level 1 heading
//   #section[Content]         -- untitled box, same colour and styling
//
// The untitled form is a continuation: it picks up the colour of the level 1
// section it sits in, i.e. the same colour as the box above it. Use it to cut
// one long box into several -- shorter boxes give the column balancer more
// places to break, so pages fill better.
#let section(first, ..rest) = {
  assert(
    rest.named().len() == 0,
    message: "#section takes no named arguments -- use #section[Title][Content] or #section[Content]",
  )
  assert(
    rest.pos().len() <= 1,
    message: "#section takes one or two content blocks -- use #section[Title][Content] or #section[Content]",
  )

  // The title bar runs edge to edge, so the inset sits on the body rather than
  // on the box around it.
  let body = if rest.pos().len() == 1 {
    heading(level: 2)[#first]
    block(width: 100%, inset: 4pt, rest.pos().first())
  } else {
    block(width: 100%, inset: 4pt, first)
  }

  // Unpainted carrier: the fill, stroke, radius and spacing all come from
  // `section_fill`, which the template hangs off the label below.
  [#block(
      width: 100%,
      breakable: true,
      body,
    )<cheatsheet-section>]
}

// Helper for section titles made of text plus a bit of inline math.
//   #section_tm[Mandelstam][$s = (sum_i p_i)^2$][ ...content... ]
#let section_tm(title, math, content) = section([#title #math], content)

#let character_table(
  classes,
  irreps,
  characters,
  class_sizes: none,
  irrep_dims: none,
  group: none,
) = {
  let n_classes = classes.len()
  let n_irreps = irreps.len()

  // Use provided irrep labels
  let irrep_labels = irreps

  // Reshape flat characters array into rows
  let char_rows = ()
  for i in range(n_irreps) {
    let row = ()
    for j in range(n_classes) {
      row.push(characters.at(i * n_classes + j))
    }
    char_rows.push(row)
  }

  // Calculate number of columns
  let n_cols = 1 + n_classes + if irrep_dims != none { 1 } else { 0 }

  // Build class sizes row (if provided) - this comes FIRST
  let size_row = ()
  if class_sizes != none {
    if irrep_dims != none {
      size_row.push([]) // empty cell for dim column
    }
    size_row.push([]) // empty cell for irrep name column
    for sz in class_sizes {
      size_row.push(text(size: 0.85em)[$#sz$])
    }
  }

  // Build header row (conjugacy class names) - this comes SECOND
  let header = ()
  if irrep_dims != none {
    header.push([]) // empty cell for dim column
  }
  // Top-left cell: group name or empty
  if group != none {
    header.push($#group$)
  } else {
    header.push([])
  }
  for cls in classes {
    header.push($#cls$)
  }

  // Build character rows
  let rows = ()
  for (i, chi_row) in char_rows.enumerate() {
    let row = ()
    if irrep_dims != none {
      row.push($#irrep_dims.at(i)$)
    }
    row.push($#irrep_labels.at(i)$)
    for chi in chi_row {
      row.push($#chi$)
    }
    rows.push(row)
  }

  // Build column alignment
  let col_align = if irrep_dims != none {
    (center, left) + (center,) * n_classes
  } else {
    (left,) + (center,) * n_classes
  }

  // Determine header row count
  let header_rows = if class_sizes != none { 2 } else { 1 }

  table(
    columns: n_cols,
    align: col_align,
    stroke: (x, y) => {
      let dim_col = if irrep_dims != none { 1 } else { 0 }
      (
        left: if x <= dim_col { none } else { 0.5pt },
        right: if x == n_cols - 1 { none } else if x < dim_col { none } else { 0.5pt },
        top: if y == 0 { none } else if y == 1 and class_sizes != none { none } else { 0.5pt },
        bottom: if y == header_rows + n_irreps - 1 { none } else if y == 0 and class_sizes != none { none } else {
          0.5pt
        },
      )
    },
    inset: 2pt,
    ..size_row,
    ..header,
    ..rows.flatten(),
  )
}


#let young_tableau(row_lengths, values) = {
  // Ensure row_lengths is an array
  let row_lengths = if type(row_lengths) == array { row_lengths } else { (row_lengths,) }

  let value_idx = 0
  let rows = ()

  // Build each row
  for row_len in row_lengths {
    let row_cells = ()

    // Add cells for this row
    for i in range(row_len) {
      if value_idx < values.len() {
        row_cells.push($#values.at(value_idx)$)
        value_idx += 1
      } else {
        row_cells.push([])
      }
    }

    rows.push(row_cells)
  }

  // Find maximum row length for table columns
  let max_cols = calc.max(..row_lengths)

  // Create table with proper alignment and separators
  let table_rows = ()

  for (idx, row_cells) in rows.enumerate() {
    let padded_row = row_cells

    // Pad row with empty cells if needed
    while padded_row.len() < max_cols {
      padded_row.push([])
    }

    table_rows += padded_row
  }

  // Determine where to place horizontal line (before last row)
  let n_rows = row_lengths.len()

  table(
    columns: max_cols,
    align: center + horizon,
    stroke: (x, y) => {
      let actual_row = y
      let actual_col = x

      // Determine if this cell should have content (within row length)
      let has_content = actual_col < row_lengths.at(actual_row)

      if not has_content {
        return none
      }

      // Draw consistent borders around every cell with content
      (
        left: 0.5pt,
        right: 0.5pt,
        top: 0.5pt,
        bottom: 0.5pt,
      )
    },
    inset: 2pt,
    ..table_rows
  )
}
