# Guaranteed assembly at a fixed stock

Fix a connected $N$-cell target $T$ and a supply of exactly $N$ loose identical cells. A rule **guarantees** $T$ if every maximal history ends in one body with that arrangement, up to rotation, translation and reflection. The rule may be designed for $T$ and $N$. A **part type** is a connected arrangement up to rotation and translation. Its reflected arrangement may have a different type. An **inventory** is a multiset of bodies, with repeated copies counted separately.

A geometric join of two types is any placement of the two bodies with disjoint interiors and at least one full-edge contact. Their union gives the output type. The [whole-body rule](rules.md) can implement choices of such joins using labels that describe each participant and its proposed contact edge.

## The catalogue theorem

**Theorem.** A set $F$ is exactly the set of part types occurring under some whole-body rule that guarantees $T$ if and only if all three conditions hold.

1. $F$ contains the singleton and at least one full-stock type. Each member fits in one hand of $T$, and each full-stock member has shape $T$.
2. Every inventory of members of $F$ with total stock $N$ and at least two bodies offers a geometric join whose output belongs to $F$.
3. Each nonsingleton member of $F$ is a geometric join of two smaller members.

In particular $T$ is achievable exactly when such a catalogue exists. A **valid factorization** is an exact catalogue as above. It is **maximal** if it has no valid strict superset.

**Proof.** If a part occurs, its entire build can be replayed using just its constituent cells. All choices depend on the participants of a join, so several of these histories can be replayed side by side, including separate copies of the same part. Every stock-$N$ inventory of occurring types can therefore be realized. It must be able to continue and eventually become $T$. Permanence embeds every part in the finish. Taking the last join that made each nonsingleton proves condition 3.

For the converse, label each exposed edge by its blank geometric part type together with the orbit of the edge under the part's proper symmetries. Use an inert label outside $F$. Two active labels are compatible exactly when placing their marked edges opposite one another is a disjoint geometric join with output in $F$. The relative motion is determined by the two edges. Choosing different representatives of the orbits gives a congruent result. Thus one symmetric table allows precisely the geometric joins staying inside $F$. Each join lowers the body count, and condition 2 forces all maximal histories to reach full stock. Condition 3 realizes all of $F$ by induction on size, with any remaining cells left unused while producing a part. ∎

**Packing test.** An inventory can be completed to $T$ by some sequence of geometric joins exactly when its members can be placed intact and disjoint to cover one hand of $T$. Necessity follows from permanence. For sufficiency, contract the placed parts in the target's whole-edge contact graph. Connectedness of that graph allows the parts to be joined in place. In a valid catalogue, *every* full-stock inventory packs jointly. Checking each part separately would miss obstructions caused by independent production.

## Pairing and powers of two

**Pairing theorem.** From an even number $m$ of identical rigid seeds $B$, a target is guaranteed exactly when some join of two copies of $B$ gives an arrangement $C$ from which the target is guaranteed with $m/2$ identical copies.

**Proof.** Choose the first allowed join of two seeds and replay it on disjoint pairs until there are $m/2$ identical copies of the resulting body. All continuations must succeed. The catalogue theorem and its proof work with a rigid seed in place of the singleton, so any decorations on the copies can be omitted when constructing a whole-body rule for $C$. Conversely, allow the chosen initial pairing, then use the catalogue rule for $C$. The stages use disjoint input types because the second never uses bodies smaller than $C$. Until all seeds have paired, there are an even number of unpaired seeds and a pair can join. All later joins can be replayed in the same order after every pair has formed, so they must finish correctly. ∎

**Corollary.** At stock $2^k$, a target is achievable exactly when it splits into two properly congruent halves, each of which has the same property recursively down to singletons. At any stock, repeatedly halving the supply reduces the question to an odd number of identical rigid seeds.

There is also a simple obstruction at every stock. Pair single cells until at most one is left. In any guaranteed target the resulting dimers must all fit, so the target's contact graph needs a matching covering all but at most one cell. At stock eight this test alone is insufficient. A seven-cell straight row with an eighth square attached above its end has such a matching. Its contact graph is a path whose only connected equal split gives a bent four and a straight four, so it cannot be guaranteed. The eight-square ring does split recursively into congruent halves. [Its maximal factorizations](../atlas/squares/stock-08.html#S08-0198)

The [four-cell square](../atlas/squares/stock-04.html#S04-0005) illustrates the difference between an exact catalogue and a maximal one. A pairing-only rule reaches just the singleton, domino and square. A rule can also allow a domino and singleton to form a bent three, and let it close with the fourth square. Both guarantee the same target. The pairing-only rule can stall on an isolated stock of three. This is the difference between the fixed-target and [all-stock](all-stocks.md) guarantees.
