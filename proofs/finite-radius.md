# Bounded-radius labels

**Theorem 1.** Fix any side count and any of the three label models. A
determinate permanent-bond rule of finite radius cannot produce the shapes
of a nonstraight zigzag window family or of a two-row ladder family.
The proof uses [the three models and the reachability lemma](rules.md)
and the [continued-family geometry](continuation.md). Shape equivalence
allows reflections.

For even side counts a blank radius-zero rule produces the straight row.
Give every singleton edge label q. On any nonsingleton, label a free edge t
if the opposite edge of its cell is bonded. Give all other edges the
inert label. Allow q–q, q–t and t–t. Inductively all bodies are rows, their
two ends are the eligible edges, and any two can join end to end.

Consequently, wherever the whole-body classification holds, the finite-radius
classification is the row for even side counts and no determinate rule for
odd side counts. The exclusion of zigzags and ladders is independent of
that classification.

**Lemma L′ (locality)**. Let X and X′ be joined into Y, and let e be a free edge of Y on a cell c of X. If no cell of X′ touches a cell of X at bond distance at most r from c in X, then N_r(Y, e) = N_r(X, e). So e has the same label in X and in Y.

*Proof.*
- Suppose a path of length at most r from c in Y uses a bond of the join. Take a shortest one, and its first join bond, from x to x′.
- Then x is within r − 1 of c in X, and x touches x′ ∈ X′. That is excluded.
- So the cells within distance r of c are the same in X and in Y, reached by the bonds of X.
- None of them touches a cell of X′, so each of their edges is in the same state in X and in Y, and their numbering is the same. ∎

### 3.1 Zigzag windows, chain case (k ≠ n/6)

Let Z be the infinite zigzag Z(n, k), with cells c_i, i ∈ ℤ. A window [a, b] is the body of cells c_a … c_b. The facts used:
- **(Z1)** Only consecutive cells touch (see [zigzag geometry](geometry.md#shape-families)). So a window's bonds are exactly its consecutive contacts, and bond distance equals index distance.
- **(Z2)** Z is carried onto itself by the translations c_i ↦ c_{i+2m} and the half-turns c_i ↦ c_{2j+1−i} (see [zigzag merges](geometry.md#why-these-families-are-merge-closed)).
- **(Z3)** Every rotation and translation g that carries a window [a, b] of at least 3 cells onto a window [a′, b′] is the restriction of a map in (Z2).
  - g preserves contacts, so by (Z1) it maps the chain onto the chain, keeping or reversing the order.
  - Turns alternate in sign along Z, and a rotation keeps the sign of a turn, while reversing the direction of travel flips it.
  - If g keeps the order, the turns at c_{a+1} and c_{a′+1} agree, so a′ − a is even, and g agrees with the translation by a′ − a on the centres of c_a and c_{a+1}.
  - If g reverses the order, a + b′ is odd, and g agrees with the half-turn that has 2j + 1 = a + b′ on those two centres.
  - A rotation and translation is fixed by where it sends two distinct points.
- **(Z4)** The mirror image of a window is a rotated copy of a window.
  - An odd window [0, 2s] is carried onto itself by a reflection ρ, because its turn word (t, −t, …, t) is a palindrome (see [zigzag geometry](geometry.md#shape-families)).
  - ρ must reverse the order, because a reflection that fixes three centres not on a line is the identity. So ρ(c_i) = c_{2s−i}.
  - For even m, take 2s = m: ρ carries [0, m − 1] onto [1, m].
- **Types.** The left end of [a, b] has type a mod 2, and its right end type (b + 1) mod 2.
  - The maps in (Z2) preserve types. So by (Z3), the type of each end of a body that is a window of at least 3 cells does not depend on how it is laid onto Z.
  - An even window has two ends of one type, which we call its *hand*. An odd window has one end of each type.
  - The *continuing edge* of an end cell is the edge shared with the next cell of Z beyond the window. It is also independent of the laying.
- **(Z5) Junctions.** Let X and Y be windows of at least 3 cells, placed so that X ∪ Y is congruent, reflections allowed, to a window. Then X and Y are the two end segments of that window. They share exactly one contact, between continuing edges at an end of X and an end of Y of equal type. The ends of X ∪ Y are the other two ends, with their types.
  - By (Z4), X ∪ Y is a rotated copy of a window [a, a + N − 1]. Lay it on Z.
  - X is connected through contacts, so by (Z1) it is a run of consecutive cells, and so is Y. So they are [a, a + p − 1] and [a + p, a + N − 1].
  - Both joined ends have type (a + p) mod 2. ∎
  - In particular, two odd windows joined at their ends of type τ make an even window of hand 1 − τ.

*Proof of Theorem 1 for these families.* Suppose a determinate rule of radius r has the windows of Z as its final shapes. By Lemma A1 and (Z4), every occurring body of at least 3 cells is a rotated copy of a window.

1. **One hand.** All occurring even windows of at least 4 cells have the same hand h.
   - If X and Y had different hands, Lemma A2(2) would give an allowed join, whose union occurs and so is a window.
   - By (Z5) the join pairs ends of equal type, but X has only ends of one type and Y only of the other.
2. **Ends of type h of odd windows never bond to each other.** Let O and O′ be occurring odd windows of at least 3 cells.
   - Lay O on Z with its type-h end on the right, and O′ as the next |O′| cells. This is possible by the maps in (Z2), and nothing overlaps.
   - Their union is an even window of hand 1 − h, and their only contact is between the two type-h continuing edges.
   - If that label pair were in the table, the join would be allowed. Carrying out both histories and then this join would make an even window of hand 1 − h occur, against step 1.
3. **Even windows built from long odd ones.** Let O₁ and O₂ occur, be odd, and have at least max(3, r + 2) cells.
   - By Lemma A2(2) some join of them is allowed, and its union E occurs, so it is a window.
   - By (Z5) and step 2, they are joined at their type-(1 − h) ends. So the ends of E are the type-h ends of O₁ and O₂.
   - Only O₁'s cell at the joint touches O₂ (Z1), and it is at bond distance |O₁| − 1 ≥ r + 1 from O₁'s free (type-h) end cell. So by Lemma L′, that end cell's continuing edge has the same label in E as in O₁. The same holds for O₂.
4. **Contradiction.** Take four such odd windows O₁ … O₄. Let E be built from O₁ and O₂, and E′ from O₃ and O₄, as in step 3.
   - By Lemma A2(2) some join of E and E′ is allowed. By (Z5) it is at an end of each, with a single contact between continuing edges.
   - Those ends are type-h ends of some O_i (i ≤ 2) and O_j (j ≥ 3), and by step 3 their labels are the labels in O_i and O_j.
   - That pair is excluded by step 2. So the join is not allowed. ∎

### 3.2 Two-row ladders (n even, 60° ≤ φ ≤ 90°)

Let Σ be the strip of the ladder construction, with cells B_i = iu (bottom row) and T_i = iu + v (top row), with columns {B_i, T_i}. Its *members* are P(2m), E(2m + 1) and T(2m + 1). [Ladder construction](geometry.md#shape-families)
- **(S0)** The contacts of Σ are row neighbours, rungs B_i–T_i, and, at 60° only, T_i–B_{i+1} (see [ladder geometry](geometry.md#why-these-families-are-merge-closed)). Every contact joins cells in the same column or in adjacent columns.
- **(S1)** Let M be a member of at least 5 cells, and g a rotation and translation with g(M) ⊂ Σ. Then g is a translation by cu with c ∈ ℤ, or a half-turn B_i ↦ T_{c−i}, T_i ↦ B_{c−i}. So g(M) is a block of full columns with at most one extra cell at an end column.
  - M has three consecutive collinear cells in a row, each touching the next.
  - By (S0), inside Σ such a triple lies in a row: a step along v, or along u − v, cannot be repeated inside Σ.
  - So g maps a row line of M onto a row line of Σ. M also has cells on its other row line, at distance sin φ. The only line of Σ's cells at that distance from a row line is the other row line. So g maps the pair of row lines onto itself.
  - If g keeps each line, it is a translation along u, by an integer since it maps cells onto cells.
  - If g swaps them, it is a half-turn about a point of the midline. Sending B_i to some T_l makes it the half-turn with c = i + l.
- **End types.** An end is *flat* if its end column is full. A step is of type *α* if its single cell is at the bottom right or the top left, as in E, and of type *β* if it is at the top right or the bottom left, as in T.
  - The maps in (S1) preserve these types. So each end of a member of at least 5 cells has a well-defined type.
  - P has two flat ends. E has one flat end and one α end. T has one flat end and one β end.
- **(S2) Junctions.** Let X and Y be members of at least 5 cells, placed so that X ∪ Y is a rotated copy of a member W. Lay it onto W ⊂ Σ. Then one of X, Y lies to the left of the other, and either:
  - *flat–flat*: the facing ends are flat and in adjacent columns; or
  - *step–step*: the facing ends are steps of the same type, sharing one column.

  All contacts between X and Y lie in at most three consecutive columns. The ends of W are the other two ends, with their types.
  - By (S1), X and Y are blocks of full columns, each with at most one extra cell at an end.
  - A column of W is full in X, full in Y, or split. When W is odd, its one-cell end column belongs to one piece. A split column holds the extra cell of each, so there is at most one. (The proof below uses only W = P.)
  - The full columns of X form an interval, and so do those of Y. So X lies on one side, and its extra cell, if any, lies at the split column or at the end of W.
  - A split column with cells B_s ∈ X and T_s ∈ Y, X on the left, makes X's right step and Y's left step both of type α. The other way round both are of type β.
  - Contacts join cells in the same or adjacent columns (S0). ∎
- **(S3)** Two odd members of one family, placed flat end to flat end, make a body with an extra cell at each end. By (S1) it is not a rotated copy of any member.
  - At 60°, for E, it is a mirror image of P.
- **The one-hand lemma is needed at every angle**, to fix one strip Σ′ below. It excludes both hands of P at oblique angles and at 60°, and mixing E and T at 90°. At 60° it is also what makes (S3) enough, since there the flat–flat union of two E bodies is a mirror image of P.

*Proof of Theorem 1 for these families.* Suppose a determinate rule of radius r has the shapes of a ladder family at angle φ as its final shapes.
- Reach is merge-closed (Theorem A(1)). So by the one-hand lemma for ladders ([continuation, Lemma 8](continuation.md#hands)), each Reach(N) is one body up to rotation and translation, and Reach or its mirror image is F_E or F_T.
- Let Σ′ be Σ or its mirror image, so that every occurring body of N cells is a rotated copy of the family's member of Σ′ with N cells. (S0)–(S3) hold for Σ′ by reflection.

1. **Flat ends of odd bodies never bond to each other.** Let O and O′ be occurring odd bodies of at least 5 cells, placed flat end to flat end: laid in Σ′ with their flat end columns adjacent, one body on each side. This includes the third contact at 60°.
   - If any contact there had its label pair in the table, the join would be allowed and the union would occur.
   - By (S3), the union is not a rotated copy of a member. So no contact there has its pair in the table.
2. **Even bodies built from long odd ones.** Let O₁ and O₂ occur, be odd, and have at least max(5, 2r + 5) cells, so at least r + 2 full columns each.
   - By Lemma A2(2) some join is allowed. Its union E occurs, so it is a rotated copy of P.
   - By (S2) the join is flat–flat or step–step. Flat–flat is excluded by step 1, since no contact there can bond. So it is step–step, and the ends of E are the flat ends of O₁ and O₂.
   - A bond moves at most one column (S0). So the cells of O₁ within bond distance r of its flat end column lie within r columns of it, and touch only cells within r + 1 columns of it.
   - Every cell of O₂ lies at least r + 2 columns from that column, so none of them touches those cells.
   - By Lemma L′, the free edges of the cells in O₁'s flat end column have the same labels in E as in O₁. The same holds for O₂.
3. **Contradiction.** Take four such odd bodies O₁ … O₄. Let E be built from O₁ and O₂, and E′ from O₃ and O₄, as in step 2.
   - By Lemma A2(2) some join of E and E′ is allowed. Its union is a rotated copy of P.
   - By (S2) it is flat–flat, since neither has a step. It joins the flat end of some O_i (i ≤ 2) to the flat end of some O_j (j ≥ 3).
   - All contacts between E and E′ lie between those two end columns (S0). They are exactly the contacts of O_i and O_j placed flat end to flat end, and by step 2 they carry the same labels.
   - By step 1 none of these pairs is in the table. So the join is not allowed. ∎
