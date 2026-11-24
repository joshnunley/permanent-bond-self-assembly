# Permanent-bond rules and their geometric reduction

Fix a side count and a stock of identical regular polygons. A state is a
collection of rigid, bond-connected bodies. A run begins with single cells.
At each step, choose two bodies and move one by rotation and translation
so that their interiors stay disjoint and a pair of their free edges
coincides. An edge is free when it has no whole-edge contact with another
cell of its body. A fixed symmetric compatibility table allows the join
if at least one pair of coinciding free edges carries compatible labels.
All compatible pairs at the interface bond. Existing bonds and cell
placements within the two bodies persist. A whole-edge contact left
unbonded stays unbonded. Labels are recomputed after each join.

A label is a function of the marked edge and the structure visible from
that edge, up to rotation and translation. The rule uses the same function
and table for every stock size. Bodies have independent placements until
they join. A join depends on its two participants and the proposed
relative placement. A run can choose any allowed join and ends when no
join is possible.

A rule is **determinate** if for every stock size N all maximal runs end
in a single N-cell body, and those bodies have the same arrangement up
to rotations, translations and reflections. Arrangement congruence
matches cells as polygons and forgets numbering, labels and bonds.
Reflections enter this comparison only. The allowed motions of a body
are rotations and translations.

## What a label sees

| Model | Numbering available to the label |
| --- | --- |
| A, whole-body numbering | Edge numbers of every visible cell |
| C, local numbering | The number of the marked edge |
| B, blank cells | Geometry and bonds, with edge numbering omitted |

Each cell starts with the same cyclic edge numbering in the numbered
models. The label is invariant under direct motions, which allows a
blank rule to distinguish handed neighbourhoods.

At radius r, the geometry includes the cells within bond-graph distance r
of the marked edge's cell, the marked edge, and the status of **every**
edge on each visible cell. The statuses are free, bonded, or in a
whole-edge contact without a bond. In particular a visible cell can show
a bond or contact to a cell beyond the radius. The numbering in the
table is then restricted to these observations. A whole-body rule,
written R∞, sees the entire body. Radius is fixed independently of
the stock size. The compatibility table can have countably many labels.

The rule classes are nested as B ⊆ C ⊆ A at either radius. Each bounded
rule is also a whole-body rule.

## Reachability lemma

Call a body *occurring* if some run, at some stock size, reaches it.
Any occurring body can be built on its own cells. Permanence makes
its assembly history a tree of joins of its own sub-bodies. Repeat
that history from its single cells. Each decision sees the same
participant bodies, so the joins are still allowed. Once the body is
complete, there is nothing left for it to join.

It follows that a rule is determinate exactly when both conditions hold.

1. Any two occurring bodies of the same size have congruent arrangements.
2. Any two occurring bodies admit an allowed join in some relative placement.

For necessity, the single-body histories force condition 1. Build the
histories of any two bodies independently in a stock of their combined
size. A join must be available when only those bodies remain. For
sufficiency, each join reduces the number of bodies by one. Condition 2
forces the last state to be a single body, and condition 1 gives its shape.

## Whole-body reduction

Let Reach(N) be the arrangements of the occurring N-cell bodies, taken
up to rotation and translation. The reachability lemma shows that
for a determinate rule, Reach is a [merge-closed family](geometry.md).
The proof applies in every model and at every radius.

Conversely, take a merge-closed family F. For each arrangement in F,
label a free edge by the proper congruence class of the entire
arrangement together with the edge's orbit under its proper
symmetries. Give every other body's edges inert labels. Two active
labels are compatible exactly when the two marked edges can be matched
with the cells on opposite sides, with disjoint interiors, and the
union's arrangement lies in F. The two edges determine the proper motion
uniquely. Choosing a different member of either edge orbit changes
the assembled arrangement by at most a proper motion, so the table is
well-defined and symmetric.

This rule reads only the blank arrangement and edge geometry, and
all contacts made by a permitted join bond. Induction from the
singleton puts each occurring arrangement in F. The merge property
then gives an allowed join between any two occurring bodies. The
reachability lemma proves determinacy. The reached arrangements
are a subset of F at each size, with equality if F has a single
proper congruence class at each size.

Thus in all three models the whole-body shape sequences are exactly
the shape sequences of merge-closed families. For each side count where
the classification is established, a whole-body rule is determinate
if and only if every pair of occurring bodies can join and its
Reach family is one of the listed families up to global reflection.
