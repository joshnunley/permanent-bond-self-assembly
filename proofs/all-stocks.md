# One guaranteed shape at every stock

A single permanent rule is **determinate** if, for every finite stock $N$, all maximal histories from $N$ loose cells finish in one body and its shape depends only on $N$. The [reachability lemma](rules.md#reachability-lemma) turns this requirement into a geometric one. If bodies of sizes $a$ and $b$ each occur, we can prepare them independently at stock $a+b$. They must be able to join. The fixed-stock problem only asks for progress when *all* bodies in the inventory total the chosen $N$.

## Classification theorem

For all odd $n\ge3$ and for even $4\le n\le22$, the possible sequences of final shapes of determinate whole-body rules on identical regular $n$-gons are exactly the following families, when geometrically admissible.

- **Straight rows**, for even $n$.
- **Alternating zigzag windows**, with bond separation $n/6\le k<n/2$.
- **Two-row ladders**, for even $n$, at available rung angles between $60^\circ$ and $90^\circ$, inclusive. At odd stock there are two choices of row for the extra cell.

This is the same in all three [label observation models](rules.md#what-a-label-sees): full numbering, the number of the marked edge only, and blank geometry. The necessity uses the richest model and the constructions use blank cells. The underlying merge-closed families have one proper congruence class at each size, with a consistent choice of hands up to reflecting the family. [Geometry and constructions](geometry.md)

For $3\le n\le23$, the [fixed-side-count proof](fixed-side-count.md) checks all prefixes through stock 16 with certified arithmetic. The [continuation proof](continuation.md) either continues a listed family or rules out an infinite continuation of a constant-turn arc. The cutoff of 16 is a search bound for these side counts. The conclusion quantifies over every stock.

For arbitrary odd $n$, the [odd-side proof](all-side-counts.md#odd-side-counts) uses contact-graph parity, a finite search independent of $n$, and the same continuation lemmas. The remaining even-side-count cases are [open](all-side-counts.md#remaining-even-side-counts).

A single settled frozen-label local rule also realizes each listed family if it may use a countable state space and settle for as long as needed. The whole-body classification limits every determinate settled local rule, since its parts can be replayed independently. [All-stock local equivalence](power-completeness.md#one-local-rule-for-every-stock)

## With a fixed observation radius

For every $n$, the [finite-radius theorem](finite-radius.md) excludes nonstraight zigzag windows and two-row ladders. A blank radius-zero rule produces rows for even $n$. Together with the classification, there is thus no finite-radius determinate rule for any odd $n$, and the row is the only finite-radius output family for even $n\le22$. Here a label is a function of its current geometric neighborhood. The separate model with [retained states and a bound on settling rounds](propagation.md) allows information from before a join.
