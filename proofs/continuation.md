# From a finite prefix to every stock size

## 1. Setting

Use the notation in [Bodies and merge-closed families](geometry.md). A cell is a regular n-gon in a definite position, so an isometry of the plane carries cells to cells and bodies to bodies. "A ≅ B" means a rotation-and-translation carries A onto B. Congruence also allows reflections. Distinct cells of a body have distinct centres, since two n-gons with the same centre overlap.

**Definition 1 (shape sequence, weak property 3).** A *shape sequence* gives, for each N ≥ 1, a congruence class S(N) of N-cell bodies, with S(1) the single cell. It has *property W at N* for the split (a, b), a + b = N, if a body X of shape S(N) is the disjoint union of a body of shape S(a) and a body of shape S(b) that touch. It is *weakly merge-closed* if it has property W at every N for every split.

If W holds for one body X of shape S(N), it holds for every body of that shape: apply the isometry.

**Lemma 1.** The shapes of a merge-closed family are weakly merge-closed.

*Proof.* Take A ∈ F(a), B ∈ F(b). Property 3 gives a member of F(a + b) that is the union of a copy of A and a copy of B that touch. Copies under ≅ are congruent. ∎

Property W forgets hands. All lock-in lemmas below assume only W, so they also apply to prefixes found by the search that hold both hands of a body somewhere.

**The split (1, M).** Property W for (M, 1) says a body of shape S(M + 1) is Y ∪ {x} with Y of shape S(M) and x one more cell.

## 2. Lock-in for chains

### The row

**Lemma 2 (row).** Let n be even and M ≥ 5. Suppose S(N) is the row of N cells for N ≤ M, and W holds at M + 1 for the splits (1, M) and (3, M − 2). Then S(M + 1) is the row of M + 1 cells.

*Proof.* Let X = Y ∪ {x} have shape S(M + 1), with Y the row c₀, …, c_{M−1}, c_i = c₀ + iu, on the line ℓ. By W for (3, M − 2), X = A ∪ B with A a row of 3 cells and B a row of M − 2 ≥ 3 cells. Let W′ be the part that contains x. It is a row of p ≥ 3 cells, and its other p − 1 ≥ 2 cells lie in Y. Its centres lie on one line, 1 apart. Two of them are on ℓ, so the line is ℓ and the centres are c₀ + tu with t running over p consecutive integers. All cells of a row are translates of one another, so x is the translate of a cell of Y by a multiple of u. Say x has centre c₀ + tu. Then t is not in {0, …, M − 1}, and a neighbour of t among the p consecutive integers is in {0, …, M − 1}. So t = −1 or t = M, and X is a row. ∎

### Constant-turn arcs

**Definition 2 (arc).** Let τ = π − 2πk/n with 1 ≤ k ≤ n − 1, k ≠ n/2. The *τ-arc* is the sequence of cells c_i (i ∈ ℤ) in which each cell touches the next and each c_i is left through the edge k steps after its entry edge. Equivalently c_{i+1} = g(c_i), where g is the rotation by τ about a point O at distance ρ = 1/(2 sin(|τ|/2)) from all the centres. Let Γ be the circle of centre O and radius ρ, and A_N = {c₀, …, c_{N−1}}. (A_N is a body only while its cells are distinct and do not overlap.)

Two facts. (i) g has finite order q, because τ is a rational multiple of π. Also c_i = c_j exactly when i ≡ j mod q, because a rotation about O fixes the centre of c₀ only if it is the identity. (ii) The reflection r₀ in the line through O and the centre of c₀ maps c_i to c_{−i}. Indeed that line bisects the angle between the two edge normals of c₀ that point to c₁ and c_{−1}. The bisector of two edge normals of a regular n-gon is a symmetry axis, so r₀(c₀) = c₀, and r₀ g r₀ = g⁻¹.

**Lemma 3 (arcs).** There is no weakly merge-closed shape sequence with S(N) = A_N for all N ≤ 7. In particular no merge-closed family has S(N) = A_N for N ≤ 7, whether or not the arc later closes into a ring.

*Proof.* Suppose there is one. We show by induction that S(N) = A_N for every N. This is absurd at N = q + 1, since A_{q+1} repeats a cell.

Let M ≥ 7 and S(N) = A_N for N ≤ M. Let X = Y ∪ {x} have shape S(M + 1) with Y = A_M. By W for (4, M − 3), the part W′ containing x is h(A_p) for an isometry h and p ≥ 4. Its other p − 1 ≥ 3 cells are cells of Y, with three distinct centres on Γ and on h(Γ). Three points lie on one circle only, so h(Γ) = Γ and h fixes O. Hence h is a rotation about O, or such a rotation after r₀. By (ii), h(A_p) = {g_α(c_{±i}) : 0 ≤ i < p} for a rotation g_α about O. One of these cells is a cell c_j of Y. Its centre shows that g_α is a power of g. So W′ = {c_s, …, c_{s+p−1}} for some s, with indices mod q.

So x = c_t, where t is not in {0, …, M − 1} mod q but t − 1 or t + 1 is. Hence t ≡ M or t ≡ −1 mod q. If M = q there is no such t, so no S(M + 1), a contradiction. Otherwise X = {c₀, …, c_M} or {c_{−1}, …, c_{M−1}} = g⁻¹{c₀, …, c_M}. Either way S(M + 1) = A_{M+1}. ∎

### Zigzag windows

Fix k with 1 ≤ k < n/2 for which Z(n, k) is defined ([shape families](geometry.md#shape-families)), with chain c_i, i ∈ ℤ, and W_N = {c₀, …, c_{N−1}}. Write a = c₁ − c₀, b = c₂ − c₁ and s = a + b. The centres are c₂ⱼ = c₀ + js and c₂ⱼ₊₁ = c₀ + a + js. The vectors a and s are nonparallel. Two centres of different parity differ by ms + a or ms + b, whose squared norms are both 1 + m(m + 1)|s|².

Call the line λ₀ through the even centres and the line λ₁ through the odd centres the *rails*. They are distinct and parallel. Cells on one rail are translates of each other by multiples of s. Two cells on different rails touch exactly when they are consecutive in the chain (1 + m(m + 1)|s|² = 1 only for m = 0, −1).

**Lemma 4 (windows).** Let M ≥ 11. Suppose S(N) = W_N for N ≤ M, and W holds at M + 1 for the splits (1, M) and (6, M − 5). Then S(M + 1) = W_{M+1}. This includes k = n/6, where the windows are not chains.

*Proof.* Let X = Y ∪ {x} with Y = W_M. Let W′ be the part of the split (6, M − 5) that contains x, with W′ = h(W_p) and p ≥ 6. The images under h of the two rails of W_p are two parallel lines, each carrying at least 3 cells of W′, equally spaced by |s|.

*W′ lies on the rails.* One of the two lines does not carry x. Its ≥ 3 cells are cells of Y with collinear centres. Three collinear centres of Y lie on one rail, because a line other than λ₀, λ₁ meets each rail at most once. So this line is a rail λ_ε. The other line carries at least 2 cells of Y, is parallel to λ_ε and different from it, so it is λ_{1−ε}. On each rail, W′ has a cell of Y, and the other cells of W′ on that rail are its translates by multiples of s (the linear part of h sends s to ±s, because the image lines are parallel to the rails and the spacing is |s|). So every cell of W′, including x, is a cell c_t of the infinite chain.

*W′ is a run of consecutive cells.* Write h(c_i) = c_{σ(i)} for 0 ≤ i < p. Since c_i touches c_{i+1} and they are on different rails, so are their images, hence σ(i + 1) = σ(i) ± 1. As σ is injective, σ(i) = σ(0) + i for all i, or σ(0) − i for all i.

So x = c_t with t outside {0, …, M − 1} and a neighbour t ± 1 inside. Hence t = M or t = −1, and X is {c₀, …, c_M} or {c_{−1}, …, c_{M−1}}. The second is congruent to the first. Let ι be the improper isometry with ι(c₀) = c₁ and ι(c₁) = c₂. It maps the chain c₀, c₁, c₂, … to the chain that starts c₁, c₂ and has the negated turn word. This gives the chain c₁, c₂, c₃, …, because the turns alternate. The same holds backwards. So ι(c_i) = c_{i+1} for all i, and ι maps {c_{−1}, …, c_{M−1}} onto {c₀, …, c_M}. ∎

## 3. Lock-in for ladders

### The strip and its forms

Let n be even and φ = 2πj/n with 60° ≤ φ ≤ 120°. Put u = (1, 0), v = (cos φ, sin φ), and let K be a cell centred at the origin with u among its edge normals. Therefore −u, ±v are edge normals too.

**Definition 3.** The *strip* Σ = Σ_φ has cells b_i = K + iu (bottom row, on the line ℓ_b) and t_i = K + iu + v (top row, on ℓ_t), i ∈ ℤ. For integers β₁ ≤ β₂ and θ₁ ≤ θ₂ the *form* [β₁, β₂; θ₁, θ₂] is {b_i : β₁ ≤ i ≤ β₂} ∪ {t_i : θ₁ ≤ i ≤ θ₂}, with *offsets* δ₁ = β₁ − θ₁ and δ₂ = β₂ − θ₂. Its size is 2(β₂ − β₁ + 1) + δ₁ − δ₂. The *E-family at φ* has the bodies

- Λ(2m) = P_m = [0, m − 1; 0, m − 1], offsets (0, 0);
- Λ(2m + 1) = E_m = [0, m; 0, m − 1], offsets (0, 1);

and Λ(1) = one cell. For 60° ≤ φ ≤ 90° this is F_E. For 90° < φ ≤ 120° it is the mirror image of F_T at 180° − φ. So *the ladder shapes are exactly the shapes Λ^φ(N), 60° ≤ φ ≤ 120°*, and no two cells of Σ_φ overlap, for all these φ (for φ > 90° it is the mirror image of the strip at 180° − φ).

A *valid form* is a form congruent to Λ(N), N its size. Sym(Σ) is the group of isometries that map the set of cells of Σ onto itself.

**Lemma 5 (symmetries and valid forms).** Put ε = 2 cos φ.

1. Sym(Σ) consists of the translations by multiples of u and the half-turns that exchange b_i and t_{c−i} (c ∈ ℤ). If ε is an integer, that is φ ∈ {60°, 90°, 120°} and ε ∈ {1, 0, −1}, it also contains the reflections x ↦ c − x in vertical lines, which map b_i ↦ b_{c−i} and t_i ↦ t_{c−i−ε}, and their products with half-turns. There are no other elements.
2. On offsets, translations act trivially, half-turns by (δ₁, δ₂) ↦ (δ₂, δ₁), vertical reflections by (δ₁, δ₂) ↦ (ε − δ₂, ε − δ₁).
3. A form is the image of some Λ(N) under a proper element of Sym(Σ) exactly when its offsets are in V⁺ = {(0, 0), (0, 1), (1, 0)}. It is the image of some Λ(N) under an improper element exactly when its offsets are in V⁻, where V⁻ = {(ε, ε), (ε − 1, ε), (ε, ε − 1)} if ε is an integer and V⁻ = ∅ otherwise. Put V = V⁺ ∪ V⁻. So, by Lemma 6.2 below, for sizes ≥ 5 the valid forms are exactly the forms with offsets in V.
4. For (δ₁, δ₂) ∈ V: |δ₁|, |δ₂| ≤ 2; |δ₁ − δ₂| ≤ 1; the size is even exactly when δ₁ = δ₂; the two rows have sizes that differ by at most 1. V⁺ and V⁻ are each invariant under (δ₁, δ₂) ↦ (δ₂, δ₁). Two forms of the same size with offsets in V are congruent, and both are congruent to Λ of that size.

*Proof.* 1. An element g of Sym(Σ) permutes the centres. A line through three centres is ℓ_b or ℓ_t, so g permutes {ℓ_b, ℓ_t}. If g keeps each line it is a translation along u or a reflection x ↦ c − x. If it swaps them it is a half-turn about a point of the mid-line or a reflection in the mid-line followed by a translation (a product of the two previous kinds). Translations must be by integers, and a reflection x ↦ c − x needs c ∈ ℤ to preserve the bottom centres. A reflection in the mid-line followed by a translation is the product of a half-turn of Sym(Σ) and a reflection x ↦ c − x, so it lies in Sym(Σ) exactly when that reflection does. The half-turn about the midpoint of b₀ and t_c sends iu to cu + v − iu, so b_i ↦ t_{c−i} and maps cells to cells because K is centrally symmetric (n even). The reflection x ↦ c − x sends the top centre (i + cos φ, sin φ) to (c − i − cos φ, sin φ) = ((c − i − ε) + cos φ, sin φ), a top centre exactly when ε is an integer. It maps cells to cells because K is symmetric in the line through its centre perpendicular to its normal u (n even).
2. Direct from these formulas. For the half-turn, [β₁, β₂; θ₁, θ₂] ↦ [c − θ₂, c − θ₁; c − β₂, c − β₁]. For the reflection, ↦ [c − β₂, c − β₁; c − ε − θ₂, c − ε − θ₁].
3. Apply 2 to the offsets (0, 0) of P_m and (0, 1) of E_m. Conversely a form with offsets (0, 0) is a translate of some P_m, one with offsets (0, 1) a translate of some E_m, and forms with the same size and offsets are translates of each other.
4. Read off from the lists: V⁻ is {(1,1), (0,1), (1,0)} at 60°, {(0,0), (−1,0), (0,−1)} at 90°, {(−1,−1), (−2,−1), (−1,−2)} at 120°. The row sizes differ by |δ₁ − δ₂|. ∎

**Lemma 6 (alignment).** Let W′ be a body congruent to Λ(p).

1. If p ≥ 6 and all cells of W′ except possibly one cell x are cells of Σ, then every cell of W′, including x, is a cell of Σ.
2. If p ≥ 5 and all cells of W′ are cells of Σ, then W′ = g(Λ(p)) for some g ∈ Sym(Σ). So W′ is a valid form with offsets in V, and any isometry h with h(Λ(p)) = W′ lies in Sym(Σ).

*Proof.* Write W′ = h(Λ(p)). The rows of Λ(p) have ⌈p/2⌉ and ⌊p/2⌋ cells. Their images, the *rows of W′*, lie on two parallel lines at distance sin φ, with centres 1 apart.

In case 1 both rows have ≥ 3 cells, and one of them does not contain x. In case 2 the longer row has ≥ 3 cells. Either way some row of W′ consists of ≥ 3 cells of Σ with collinear centres. Three collinear centres of Σ lie in one row of Σ, because any other line meets ℓ_b and ℓ_t at most once each. So that row of W′ lies on ℓ_b or ℓ_t. The other row of W′ contains a cell of Σ (at least 2 in case 1, all of them in case 2), lies on a parallel line at distance sin φ, hence on the other of ℓ_b, ℓ_t.

So h maps {ℓ_b, ℓ_t} to itself. On each of the two lines, h maps the centres of Σ (a set of points 1 apart) to a set of points 1 apart that contains a centre of Σ, hence onto the centres of Σ on the image line. Finally all cells of Σ are translates of K, and h maps at least one cell of Σ to a cell of Σ, so the linear part of h preserves K, and h maps every cell of Σ to the cell of Σ with the image centre. Hence h ∈ Sym(Σ), and W′ = h(Λ(p)) ⊂ Σ. ∎

### Lock-in

**Lemma 7 (ladders).** Let M ≥ 11. Suppose S(N) = Λ^φ(N) for all N ≤ M, and W holds at M + 1 for the splits (1, M) and (6, M − 5). Then S(M + 1) = Λ^φ(M + 1).

*Proof.* Let X = Y ∪ {x} have shape S(M + 1), with Y = Λ(M) ⊂ Σ. By W, X = A ∪ B with A of shape Λ(6) and B of shape Λ(M − 5), M − 5 ≥ 6. Each of A, B has all its cells in Σ except possibly x. By Lemma 6 (part 1, then part 2), x is a cell of Σ and A, B are valid forms. By Lemma 5.4 each row of A and of B has at least 3 cells.

Let L be the one of A, B whose bottom row is on the left, and R the other: L = [l₁, l₂; l₁′, l₂′], R = [r₁, r₂; r₁′, r₂′] with l₂ < r₁.

*No crossing: l₂′ < r₁′.* Otherwise r₂′ < l₁′. Write |L_b|, |R_b|, |L_t|, |R_t| for the numbers of bottom and top cells. From δ₁(L) ≥ −2, l₁′ ≤ l₁ + 2. From δ₂(R) ≤ 2, r₂ ≤ r₂′ + 2 ≤ l₁′ + 1 ≤ l₁ + 3. But r₂ = r₁ + |R_b| − 1 ≥ l₂ + |R_b| = l₁ + |L_b| + |R_b| − 1. So |L_b| + |R_b| ≤ 4. In the same way, from δ₁(R) ≤ 2 and δ₂(L) ≥ −2: l₂′ ≤ l₂ + 2 ≤ r₁ + 1 ≤ r₁′ + 3, while l₂′ = l₁′ + |L_t| − 1 ≥ r₂′ + |L_t| = r₁′ + |R_t| + |L_t| − 1. So |L_t| + |R_t| ≤ 4. Then X has at most 8 cells, but it has M + 1 ≥ 12. (This step uses only |δ| ≤ 2 and that L ∪ R has more than 8 cells.)

*No gaps: r₁ = l₂ + 1 and r₁′ = l₂′ + 1.* The rows of Y are intervals, and X has one more cell. If the bottom indices of X did not form an interval, they would be the bottom interval of Y and the single index of x, and one of L, R would have a bottom row of one cell. Rows have at least 3 cells. The same holds for the top row.

So X is the form [l₁, r₂; l₁′, r₂′], and with δ := l₂ − l₂′ we have δ₂(L) = δ = δ₁(R), δ₁(X) = δ₁(L), δ₂(X) = δ₂(R).

One of L, R is A, of even size 6, so its two offsets are equal (Lemma 5.4). If L = A, then δ₁(X) = δ₁(L) = δ, so the offsets of X are (δ₁(R), δ₂(R)) ∈ V. If R = A, then δ₂(X) = δ, so the offsets of X are (δ₁(L), δ₂(L)) ∈ V. By Lemma 5.4, X is congruent to Λ(M + 1). ∎

The E-family at 60° is the mirror image of Z(n, n/6), so Lemmas 4 and 7 overlap there.

### Hands

**Lemma 8 (ladders: one hand, and the same hand throughout).** Let F be merge-closed with S(N) = Λ^φ(N) for every N. Then F(N) holds one body for every N, and F or F̄ is the E-family at φ. F or F̄ is F_E at φ if φ ≤ 90°, and F_T at 180° − φ if φ > 90°.

*Proof.* *Splits of a ladder inside the strip.* Let a, b ≥ 5 and let Λ(a + b) ⊂ Σ be the union of disjoint bodies C, D congruent to Λ(a), Λ(b). By Lemma 6.2 both are valid forms. Let L be the one whose bottom row is on the left and R the other. They do not cross, by the step in the proof of Lemma 7: a crossing pair has at most 8 cells, and a + b ≥ 10. They leave no gap: the rows of Λ(a + b) are intervals, and each row of L and of R is a non-empty interval. So the offsets are

  L: (δ₁(Λ), δ), R: (δ, δ₂(Λ)), where (δ₁(Λ), δ₂(Λ)) = (0, 0) if a + b is even and (0, 1) if a + b is odd. (∗)

*One hand.* Suppose F(N) = {A, Ā} with A ≇ Ā. Then N ≥ 5. Λ(1) and Λ(2) have no hand. Λ(3) is symmetric in the bisector of the two contact directions at its middle cell (or is a triangle). Λ(4) is a rhombus, symmetric in a diagonal. Property 3 for the pair (A, Ā) gives a member T of F(2N) that is the union of a copy C of A and a copy D of Ā. Move T by an isometry κ onto P_N ⊂ Σ. Then κ(C), κ(D) are congruent but not ≅. By (∗) they are forms with offsets (0, δ) and (δ, 0) and the same size. A half-turn in Sym(Σ) swaps the offsets, and forms with equal size and offsets are translates. So κ(C) ≅ κ(D). Contradiction.

*The same hand.* Now F(N) = {A_N}. For each N choose an isometry κ_N with κ_N(A_N) = Λ(N), and let σ_N = ±1 be its sign (+1 if proper). Call N *handed* if Λ(N) has a hand. Then σ_N does not depend on the choice. For a valid form f of size N let τ(f) be the set of signs of isometries carrying Λ(N) onto f. By Lemmas 6.2 and 5.3, +1 ∈ τ(f) iff the offsets of f are in V⁺, and −1 ∈ τ(f) iff they are in V⁻. N is handed iff τ(f) is a single sign for one (hence every) valid form f of size N. (All this is for N ≥ 5, where Lemma 6.2 applies. We use a, b ≥ 5 below.)

Let a, b ≥ 5. Property 3 gives A_{a+b} = C ∪ D with C = g₁(A_a), D = g₂(A_b), g₁, g₂ proper. Apply κ = κ_{a+b}. Then κ(C) = κ g₁ κ_a⁻¹ (Λ(a)), so σ_{a+b} σ_a ∈ τ(κ(C)), and likewise σ_{a+b} σ_b ∈ τ(κ(D)); and κ(C), κ(D) are the L and R of (∗) in some order.

- a + b even. L and R have offsets (0, δ) and (δ, 0). V⁺ and V⁻ are invariant under swapping the offsets, so τ(L) = τ(R). If a and b are handed, this set is a single sign, and σ_a = σ_b.
- a + b odd, and both a and b handed. Sizes of even parity are handed exactly when (0, 0) ∉ V⁻, or ε ≠ 0. Odd sizes are handed exactly when (0, 1) ∉ V⁻, that is ε ≠ 1. So here ε is not an integer or ε = −1, and V⁺ ∩ V⁻ = ∅. The part of even size has equal offsets (Lemma 5.4). If it is L, then δ = 0, L has offsets (0, 0) and R has (0, 1), both in V⁺ and outside V⁻. Thus σ_{a+b} σ_a = +1 and σ_{a+b} σ_b = +1. If it is R, then δ = 1 and R has offsets (1, 1), which is not in V for these ε. So σ_a = σ_b.

Hence σ_N takes one value σ on all handed N ≥ 5. Handed sizes of one parity are linked by the first case (a, b ≥ 5 of that parity), and the two parities by the second. There are no handed N < 5. If σ = +1, each A_N ≅ Λ(N), so F is the E-family. If σ = −1, F̄ is. ∎


### Isolation

**Theorem 9 (isolation).** Let 𝓛 be a listed family for n-gons: the row (n even), a window family Z(n, k), or the E-family at some φ = 2πj/n in [60°, 120°] (n even). Let F be merge-closed, and suppose S(N) is the shape of 𝓛(N) for every N ≤ 11. Then S(N) is the shape of 𝓛(N) for every N, F(N) holds one body for every N, and F or F̄ is 𝓛.

*Proof.* Lemma 1 gives W. Induct with Lemma 2, 4 or 7 to obtain the listed shapes at every size. For ladders and for Z(n, k) at k = n/6, Lemma 8 fixes the hands. A straight row is achiral.

For a zigzag with k ≠ n/6, the contact graph is a path. If a size N admits both hands H and H̄, the pairs (H, H), (H, H̄) and (H̄, H̄) must all merge at size 2N. Each body at size 2N supplies only one pair, via its single middle cut, and at most two bodies are available. Thus each size has one hand. Odd windows and the dimer are achiral. At even size N ≥ 4, a window's alternating turn word has even length. Reversing its order and negating all turns preserves its initial sign, so this sign specifies its hand. For two even sizes at least four, the cut producing those sizes occurs at an even index. Both halves of the cut then have the same initial turn sign. This forces agreement at every even size, which gives the zigzag family or its global reflection. ∎
