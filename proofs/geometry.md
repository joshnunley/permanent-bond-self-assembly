# Bodies and merge-closed families

## Definitions

- **Cells.** Fix n ≥ 3. A *cell* is a closed regular n-gon with inradius (apothem) 1/2. Cells are placed in the plane by rotations and translations only.
- **Touching.** Two cells *touch* if a whole edge of one coincides with a whole edge of the other. Their centres are then exactly 1 apart, and the line through the centres is perpendicular to the shared edge. Corner contacts and partial edge contacts do not count. Cells need not lie on any lattice.
  - Conversely, two cells with disjoint interiors whose centres are exactly 1 apart touch. Their incircles meet at one point p. If p were inside either cell, the cells' interiors would meet near p. So p is on the boundary of each cell and on its incircle, which makes it the midpoint of an edge of each. Both edges lie on the common tangent at p and have the same length, so they coincide. The code finds contacts this way, by centre distance.
- **Bodies.** A *body* is a finite set of cells with pairwise disjoint interiors, whose *contact graph* is connected. The contact graph has the cells as vertices, and an edge between two cells when they touch.
- **Chains.** A body is a *chain* if its contact graph is a path: its cells can be ordered so that each touches only the cells just before and after it.
- **Comparing bodies.**
  - A ≅ B means some rotation-and-translation carries the cells of A onto the cells of B.
  - A and B are *congruent* if some isometry (reflections allowed) does.
  - The *mirror image* Ā of A is its image under a reflection in a line. It is well defined up to ≅.
  - A body has a *hand* if it is not ≅ to its mirror image.
  - Isometries preserve touching, so congruent bodies have isomorphic contact graphs.
- **Merge-closed family** (the central definition). A family F gives, for each N ≥ 1, a nonempty set F(N) of N-cell bodies, taken up to ≅, such that:
  1. F(1) consists of the single cell.
  2. the bodies in each F(N) are congruent to one another. So F(N) holds one body, or two bodies that are mirror images of each other (a body with a hand, and its mirror).
  3. for all a, b ≥ 1 and **all** A ∈ F(a), B ∈ F(b): B can be placed against A, by rotation and translation only, so that:
     - interiors stay disjoint.
     - some cell of B touches some cell of A.
     - the union is ≅ to a member of F(a + b).

  Equivalently, for 3: for all A ∈ F(a), B ∈ F(b), some member of F(a + b) is the disjoint union of a copy of A and a copy of B (copies under ≅) that touch somewhere.
- **Notation.**
  - S(N) is the common shape of the bodies in F(N): their class under congruence. Statements like "S(N) is a chain" do not depend on which member is meant.
  - "F(N) holds one body" means |F(N)| = 1 as a set taken up to ≅.
  - The *mirror image* of a family F is the family F̄ with F̄(N) = {B̄ : B ∈ F(N)}. It is merge-closed whenever F is.
## Shape families

Each family below is merge-closed for every stock size.
- **The row** (n even). N cells in a straight line, each touching its neighbours through a pair of opposite edges.
- **Zigzag windows.**
  - Take k with 1 ≤ k < n/2. Let c₀, c₁, c₂, … be the infinite chain in which each cell touches the next, and each inner cell's two contact edges are k edges apart.
  - Going round the cell anticlockwise, the exit edge is k steps before the entry edge at c₁, k steps after it at c₂, and so on alternately. So the turns alternate.
  - The *window family* Z(n, k) has F(N) = {c₀ … c_{N−1}}.
  - The admissible parameters are exactly n/6 ≤ k < n/2. The three-cell calculation below gives the lower bound. If c₀, c₁, c₂ have disjoint interiors, then the whole chain does.
    - Write a = c₁ − c₀, b = c₂ − c₁ and s = a + b = c₂ − c₀. Because the turns alternate, c₂ⱼ = c₀ + js and c₂ⱼ₊₁ = c₀ + js + a.
    - So the difference of two centres is ms when the indices have the same parity, and ms + a or ms + b otherwise, for an integer m. (For example, c₂ⱼ − c₂ᵢ₊₁ = (j − i − 1)s + b.)
    - Let ψ be the angle between a and b. Then |s|² = 2 + 2 cos ψ and s·a = s·b = 1 + cos ψ = |s|²/2. Hence |ms + a|² = |ms + b|² = 1 + m(m + 1)|s|².
    - Each three consecutive cells cᵢ, cᵢ₊₁, cᵢ₊₂ are a copy of c₀, c₁, c₂ or of its mirror image, so cells two apart do not overlap. Their centres are |s| apart, and their incircles (radius 1/2) are disjoint, so |s| ≥ 1.
    - Every other pair of cells that are not neighbours has centres at least √3 apart: m|s| ≥ 2 when |m| ≥ 2, and 1 + m(m + 1)|s|² ≥ 3 when m(m + 1) ≥ 2.
    - For n ≥ 4 the circumradius R of a cell is at most 1/√2, so 2R < √3, and such cells are disjoint. For n = 3 the only window family is Z(3, 1), a strip of the triangular tiling.
  - For k = n/6, cells two apart also touch, so the bodies are not chains.
  - The chain c and its mirror image give mirror-image families.
  - Handedness of windows.
    - A window with an odd number of cells has no hand. Write its turns at c₁, …, c_{N−2} as a sequence of signs. With k fixed, this sequence determines the window up to ≅. Reflecting the window negates every turn. Reading it backwards reverses the sequence and also negates every turn. So the mirror image, read backwards, has the reversed sequence. For an odd window the alternating sequence has odd length, so it is its own reverse, and the mirror image is the same body.
    - So F(N) holds one body at every odd N, whichever hand of the chain is used.
    - Windows of 2 cells have no hand either.
    - A window with an even number N ≥ 6 of cells has a hand. Its centres form a 2 × (N/2) grid, spanned by s and a. The grid's outline is a parallelogram with sides (N/2 − 1)|s| ≥ 2 and 1, so it is not a rhombus. Since a·s = 1 + cos ψ > 0, its angle is not 90°, so it is not a rectangle. It is not degenerate: a is parallel to s only if ψ = 0 (exit edge opposite the entry edge, k = n/2, excluded) or ψ = π (c₂ = c₀, impossible). If the window were ≅ to its mirror image, some reflection-type isometry would carry it onto itself. Such an isometry fixes the centroid of the centres, so it is a reflection in a line, not a glide. It preserves the outline, and a parallelogram that is neither a rhombus nor a rectangle has no mirror symmetry.
    - So every window family has a hand at some N.
- **Two-row ladders** (n even).
  - Orient the cell so that its edge normals point in the directions 2πe/n (e = 0, …, n − 1). Let u = (1, 0), the normal of edge 0. Since n is even, −u is also a normal. Let v = (cos φ, sin φ) be another edge normal, φ = 2πj/n, with 60° ≤ φ ≤ 90°.
  - The *strip* has cells at iu (bottom row) and iu + v (top row), i ∈ ℤ, all translates of this one cell. Cells in a row touch, and so do rung partners iu and iu + v.
  - Its bodies:
    - P(2m) = bottom 0 … m − 1 together with top 0 … m − 1.
    - E(2m + 1) = P(2m) plus bottom cell m.
    - T(2m + 1) = P(2m) plus top cell m.
  - F_E has F(N) = {P(N)} at even N and {E(N)} at odd N. F_T is the same with T. Both are proved merge-closed for every N and every such φ.
  - Notes on the angle:
    - Below 60° the rows overlap.
    - Angles above 90° give nothing new. Reflecting in the vertical axis turns the strip at φ into the strip at 180° − φ. It turns E(N) at φ into T(N) at 180° − φ. The reflection moves the extra bottom cell to the left end, and a half-turn of the strip (which swaps the rows and reverses them) puts it at the right end of the top row. So F_T at 180° − φ is the mirror image of F_E at φ, and vice versa.
    - At 90° (n divisible by 4), E and T are mirror images, so F_T is the mirror of F_E.
    - At 60° (n divisible by 6), top i also touches bottom i + 1. Read as the chain bottom 0, top 0, bottom 1, top 1, …, the strip is the mirror image of the chain of Z(n, n/6), with the turn convention above. So F_E at 60° is the mirror image of Z(n, n/6).
    - F_T at 60° is a different family: T(3) is a chain of three cells, while E(3) is a triangle of three mutually touching cells.
    - Strictly between 60° and 90°, E and T are not congruent.

## Admissible zigzags

Put ε = π/n and orient the first step of a zigzag along the outward
normal of angle zero. The next step makes angle π − 2kε with it.
Cells 0 and 2 are translates of one another, since each step is a
half-turn about the common edge midpoint. Their displacement has
length d = 2 sin(kε) and angle (n/2 − k)ε. The incircles overlap
when d < 1, so k ≥ n/6 is necessary.

For even n, the distance between opposite supporting lines of a
cell in a direction normal to an edge is 1. At k = n/6, the
displacement has length 1 along the edge normal at angle π/3.
For k > n/6 we have 6k − n ≥ 2. Some edge normal makes angle
at most ε with the displacement. Its component in that direction
is at least
$2\sin(\pi/6+\varepsilon/3)\cos\varepsilon$.
This lower bound is at least 1. It equals
$\sin(\pi/6+4\varepsilon/3)+\sin(\pi/6-2\varepsilon/3)$,
a concave function with value 1 at both endpoints of [0, π/4].
The supporting lines then separate the interiors.

For odd n ≥ 5, the edge normals together with their opposites are
spaced by ε. The displacement direction is halfway between two of
them. In either of those directions the width of a cell is
$1/2+1/(2\cos\varepsilon)$, with an edge on one side and a vertex
on the other. The component of the displacement along the direction
is $d\cos(\varepsilon/2)$. The required inequality that this
component be at least the width is equivalent to
$d\cos\varepsilon \ge \cos(\varepsilon/2)$.
Use 6k − n ≥ 1 and the monotonicity of sine on [0, π/2].
The difference
$2\sin(\pi/6+\varepsilon/6)\cos\varepsilon-\cos(\varepsilon/2)$
vanishes at ε = 0 and ε = π/5. It is concave between them.
Indeed, the sine product is the sum of
$\sin(\pi/6+7\varepsilon/6)$ and $\sin(\pi/6-5\varepsilon/6)$.
In the second derivative the first contributes at most −49/72, the
second contributes at most zero, and the cosine contributes at most 1/4.
The difference is therefore nonnegative throughout. For n = 3 the
sole choice k = 1 gives a strip of the triangular tiling. ∎

## Why these families are merge-closed

The straight row admits every split by a cut between two cells.

For a zigzag, a half-turn about the midpoint of the bond between cells
$j$ and $j+1$ maps cell $i$ to cell $2j+1-i$. The half-turn carries
the two cells at the bond onto one another, then carries the remaining
cells onto the alternating continuation. Shifting the window by an even
number of indices is a translation, from the product of two such half-turns.
In the window from index 0 through a + b − 1, put an even-sized part first
whenever there is one. The second part then starts at an even index and
is a translate of the prescribed window. If both sizes are odd, the
half-turn about the midpoint of the bond at index (a + b − 2)/2 maps
the suffix of size b onto the window from index 0 through b − 1. In each case the two parts touch
at the cut. This proves the merge property with a consistent proper
congruence class at every size.

For ladders, the following argument uses the circumradius
R = 1/(2 cos(π/n)).

- **No overlap.**
  - Two translated copies K and K + w of a regular n-gon K with apothem 1/2 overlap exactly when ⟨w, ν⟩ < 1 for every edge normal ν. This holds because K + w meets the interior of K exactly when w lies in the interior of K − K. For even n, K is centrally symmetric, so K − K = 2K, the n-gon with apothem 1. (This argument uses even n.)
  - Cells in one row are 1 apart along an edge normal, so they share an edge. So do rung partners iu and iu + v.
  - The closest other pair is top cell i and bottom cell i + 1.
    - Here w = u − v. Its length is 2 sin(φ/2), and it points at angle φ/2 − 90° = (j − n/2)·π/n.
    - If j − n/2 is even, that direction is an edge normal. Otherwise it is π/n from one.
    - So the pair does not overlap if 2 sin(πj/n)·cos(π/n) ≥ 1, or just 2 sin(πj/n) ≥ 1 in the even case.
  - At 60° (n = 6j), j − n/2 = −2j is even and the distance is exactly 1, so the cells share an edge. These are the additional contacts in Z(n, n/6).
  - Above 60°:
    - 6j − n is even and positive, so 6j − n ≥ 2 and πj/n ≥ π/6 + π/(3n).
    - Put ε = π/n and f(ε) = 2 sin(π/6 + ε/3)·cos ε = sin(π/6 + 4ε/3) + sin(π/6 − 2ε/3).
    - On [0, π/4] both arguments lie in [0, π/2], so f is concave there. Also f(0) = f(π/4) = 1, so f ≥ 1 on [0, π/4].
    - Since n ≥ 4, ε ≤ π/4, so the pair does not overlap. Moreover f > 1 unless n = 4. (This uses sin(πj/n) ≥ sin(π/6 + ε/3), which holds because πj/n = φ/2 ≤ π/4, where sine is increasing.)
  - Every other pair is farther apart. 2R = 1/cos(π/n) ≤ √2 for n ≥ 4, and:
    - bottom i and top i + k (k ≥ 1) are √(k² + 2k cos φ + 1) ≥ √2 apart.
    - top i and bottom i + k (k ≥ 2) are at least √3 apart.
    - cells two or more apart in one row are at least 2 apart.
  - Below 60°, top i and bottom i + 1 are 2 sin(φ/2) < 1 apart, so they overlap.
- **Merges.** The strip is carried onto itself by the translations iu, and by the half-turn about the midpoint of a rung. The half-turn swaps the rows and reverses them: it maps bottom cell i to top cell −i.
  - a and b even: P(a) followed by P(b) is P(a + b).
  - a odd, b even: P(b) followed by E(a) (or T(a)) is E(a + b) (or T(a + b)).
  - a and b odd, a = 2m + 1, b = 2q + 1:
    - E(a) is P(2m) plus bottom cell m.
    - The half-turn of E(b) is P(2q) plus a top cell at its left end.
    - Place it with its parallelogram in columns m + 1 … m + q and its extra cell at top m. Together the two fill P(a + b).
    - For T, swap the rows: T(a) has its extra cell at top m, and the half-turned T(b) puts its extra cell at bottom m.
  - In each case the two parts share an edge (for example, bottom m and bottom m + 1).
- **E and T differ.** E(3) is a bent 3-chain with angle φ at its middle cell, and T(3) one with angle 180° − φ. They are congruent only at 90°. (At 60°, E(3) also closes into a triangle of three cells.) ∎
