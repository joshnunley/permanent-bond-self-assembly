# Exhaustive classification for 3 ≤ n ≤ 23

**Theorem.** For each integer n in this range, a merge-closed family of
regular n-gon arrangements has one proper congruence class at each stock
size. Up to reflecting the entire family, it is a straight row, an
admissible zigzag window family, or a two-row ladder family, with the
parameters in [the geometric definitions](geometry.md).

The proof combines a finite exhaustive search with
[continuation and handedness lemmas](continuation.md). The stock cutoff
is 16 for the computation. The theorem quantifies over every stock size.

## Arithmetic and completeness

Normalize the centre of a cell to 0 and one outward edge normal to 1.
Put ζ = exp(πi/n). All centres in a contact-connected arrangement then lie
in the integer ring ℤ[ζ]. Across an edge contact the distance is 1 and
the normal reverses. Therefore every step in a contact path is a power
of ζ, and the unoriented polygon at a centre with a contact is uniquely
oriented by that normal. For congruence of an arrangement with
at least two cells, its centre set suffices. A singleton is unique.

A candidate only needs to pass the weaker overlap condition that centres
are at least 1 apart. Every genuine polygon arrangement passes it because
the polygons contain open incircles of radius 1/2. At distance 1,
disjoint genuine polygons necessarily touch across an entire edge, as in
[the geometric definitions](geometry.md#definitions). In the relaxed
search, any contact must also be along an edge normal of both cells.
Working in this larger class can retain extra prefixes while preserving
all genuine prefixes.

Represent a centre by its integer coefficients in the power basis
modulo the cyclotomic polynomial Φ₂ₙ. Coincidence is equality of integer
vectors. A difference of modulus 1 in ℤ[ζ] is a power of ζ. Indeed,
its product with its conjugate is 1. The cyclotomic Galois group is
abelian, so all its conjugates also have modulus 1. Kronecker's theorem
makes it a root of unity, and the roots of unity in a cyclotomic field
of even index 2n are the 2n-th roots. For other differences, enclose the
real and imaginary components with outward integer bounds obtained from
interval sine and cosine. Compare the resulting enclosure for squared
length against 1. If that interval does not separate them, increase the
working precision in interval arithmetic. The search raises an error if
a sign cannot be certified. The n-gon search also checks the integer
coordinate budget before using int64 operations.

To compare proper congruence classes, map each ordered contact to the
ordered pair (0, 1), transform the remaining centres by the same direct
motion, sort them and take the lexicographically least tuple. Any direct
congruence matches the ordered contacts and thus the tuples. To test
whether a candidate is a union of two prescribed arrangements, try the
images of an ordered contact of one part at each ordered contact of the
candidate, then compare the complement with the other part. The test
handles singleton parts by placing the other part, or checking the dimer.

## Search and continuation

Start with the singleton. Every size-N arrangement in a merge-closed
family can be obtained by adding one cell to an arrangement from size
N − 1. Conversely, every prior arrangement must admit at least one such
extension. Generate extensions by trying free edges. For each new shape,
try the two proper hands separately and the two together, or one choice
if the body is achiral. Retain a choice only if every split a + b = N and
every prescribed pair of hands at a and b can be realized in the size-N
choice. By induction, the prefix of any merge-closed family occurs in
this search.

Compare the shape sequence of each surviving prefix through size 16,
including all its earlier shapes, with the row, zigzag windows, ladders
and constant-turn arcs. The checker reports zero unmatched
prefixes for the stated side counts.
[Exact checker](../code/exact_n.py) · [Recorded runs](../results/)

A shape sequence that starts as a constant-turn arc cannot continue
to all stocks, even if an initial segment closes into a ring. This is
Lemma 3 of [the continuation proof](continuation.md#constant-turn-arcs).
Every other surviving shape sequence matches one of the families
through size 11. Theorem 9 of [the continuation proof](continuation.md#isolation)
forces that family at every stock, and fixes its proper congruence
classes up to global reflection. The constructions in
[the geometry proof](geometry.md#why-these-families-are-merge-closed)
give the converse. ∎
