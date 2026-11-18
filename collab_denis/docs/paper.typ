#let title = [bay stat thema]
#let authors = ("Nicolas Volken", "Denis Titov", "Ian Röthlisberger", "Sophia Kugler")
#let course = [Bayesian Statistics, HS2026]

#set document(title: title, author: authors)
#set text(font: "New Computer Modern", size: 11pt, lang: "en")
#set par(justify: true)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")

// Placeholder until the image files are added to the repo
#let todo-figure(caption, height: 6cm) = figure(
  rect(width: 60%, height: height, fill: luma(230), stroke: none)[
    #align(center + horizon, text(fill: luma(120))[image missing])
  ],
  caption: caption,
)

// ---------------------------------------------------------------------------
// Title page
// ---------------------------------------------------------------------------

#page(header: none, footer: none)[
  #v(4cm)
  #align(center)[
    #text(size: 17pt, title)
    #v(1.5em)
    #text(size: 12pt)[#authors.join(", ") \ #course]
    #v(2.5em)
    #text(weight: "bold")[Abstract]
  ]
  #pad(x: 1cm)[
    blabal

    notes: → selection affects results (biased measurement?) → try different eccentricities, compare with other data from similar planet/star systems
  ]
]

// ---------------------------------------------------------------------------
// Running header/footer for the rest of the document
// ---------------------------------------------------------------------------

#set page(
  header: context {
    let page-no = counter(page).display()
    // Like the LaTeX version: only unnumbered pages (Contents, References) get a label
    let on-page = query(heading.where(level: 1)).filter(h => h.location().page() == here().page())
    let label = if on-page.len() > 0 and on-page.first().numbering == none {
      upper(on-page.first().body)
    } else { [] }
    grid(columns: (1fr, auto), label, page-no)
    v(-0.6em)
    line(length: 100%, stroke: 0.4pt)
  },
  footer: grid(columns: (1fr, auto), [ETH Zürich], [Sophia Kugler]),
)

#outline()
#pagebreak()

#image("/assets/image.png")

= Introduction
blabla

= Theoretical Background

== Exoplanet Detection
There are many ways to detect exoplanets. The most used method however is the transit method. When a planet moves in front of a star, the intensity gets reduced by a factor proportional to:
$ Delta I tilde R_P^2 / R_S^2 $
where $R_P$ is the radius of the planet and $R_S$ the radius of the star. The transit probability is \*\*\*; so the probability that we observe a transit is less likely for smaller planets further out. Thus the transit method favours large, close-in planets.

== Host Star
copied from wikipedia: Host Star (HD 137010)

- Star Type: K-type main-sequence dwarf (K3.5V)
- Mass Radius: 73
- Temperature: 4,770 K (about 1,000 K cooler than the Sun)
- Stellar Flux: Receives roughly 29

#link("https://exoplanet.eu/catalog/hd_137010_b--11844/")

== MCMC

== more bay stat theory as this is not a project on astrophysics but on bayesian statistics

#pagebreak()

#image("/assets/image-1.png", width: 60%)

= Methods
blabla

= Results and Discussion
blabla

= Conclusion and Outlook
blabla

#pagebreak()

#heading(numbering: none)[References]
