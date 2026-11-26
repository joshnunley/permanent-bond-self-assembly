# Power-completeness of frozen-label local rules

## The rule

All cells start in the same state and have the same cyclic port numbering. A fixed symmetric table specifies which labels can bind. A cell's state supplies the labels on its free ports. When a join forms a bond, the labels actually used at its ends are stored at those ports for the lifetime of the bond.

A cell updates by reading its own immutable bond records and, for each bonded port, the neighboring cell's state, the number of its facing port, and the frozen label at that port. The only mutable inputs are the neighbors' states. The facing label can equivalently be copied into the neighbor's state. Communication uses bonds. An unbonded contact does not exchange messages.

After each join, local updates settle to a fixed point before another join. A step can update any nonempty subset of cells, reading states from before that step. A schedule is fair if a cell whose update would change its state is eventually selected while that change remains enabled. Every reachable settling outcome counts. A fair history that fails to settle fails the guarantee. Cells hold their states from one join to the next.

## The theorem

**Theorem.** Fix the polygon and stock $N$. Allow the local rule, its state space and its settling time to depend on $N$. With a general symmetric binding table, frozen-label local rules realize exactly the same guaranteed targets and exact geometric catalogues as full-shape whole-body rules. We call this **power-completeness**.

More precisely, any whole-body rule whose edge labels are a function of the current numbered geometry, bond graph and frozen bond labels has a frozen-label implementation with the *same binding table, settled free labels and allowed joins* at that stock.

**Proof of the implementation.** Work in each cell's numbered frame. The two facing port numbers determine the relative placement across a bond. At a cell $v$, let $V_0(v)$ be its own cell and the records of its incident bonds, including the facing ports and frozen labels. At each update, compute that record and then the levels

$$
V_j(v)=V_0(v)\;\cup\!
\bigcup_{w\text{ bonded to }v}
P_{v\leftarrow w}V_{j-1}(w),
\qquad 1\le j\le N-1 ,
$$

where $P_{v\leftarrow w}$ converts the neighbor's coordinates to $v$'s and the neighbor's reports come from the start of the update step. Use the top level to compute the desired free labels, with an inert fallback if its report is incomplete or inconsistent.

After a join, level zero has fixed input. Once level $j-1$ is correct everywhere, any update corrects level $j$ and it stays correct until the next join. Fairness therefore stabilizes the levels in order. Each fixed level contains the reports from within that many bond steps. The top level covers the entire body, since it is connected and has at most $N$ cells. From the full placements one can also recognize unbonded contacts. The state space can be finite because the reports have bounded depth, only finitely many relative placements per bond step, and finitely many reachable frozen labels at fixed stock. Initialize the singleton with its correct report and labels. At each settled join decision both rules show the same labels, so they permit the same placements and create the same bonds. Induction proves the implementation claim. ∎

**Proof of the converse.** Choose a settled witness history for every part type that occurs under a successful local rule. The witnesses can be replayed independently with their states, numbering and update schedules. A settled body can wait while another is built. Thus every stock-$N$ inventory from the geometric catalogue can be realized in some concrete states. Its continuation under the successful local rule gives a geometric join within the catalogue. A last join gives the production condition for each nonsingleton. The [catalogue theorem](fixed-target.md#the-catalogue-theorem) now gives a blank whole-body rule with the same exact set of part types. ∎

## One local rule for every stock

**Proposition.** If the common local rule may have a countable state space and can settle for as long as needed between joins, it has the same determinate shape families as the whole-body model.

**Proof.** To simulate a whole-body rule without choosing a stock bound, let each cell hold a single finite report. Store placed cells and positive bond records in the report. At every update take the union of the cell's own freshly read record and all neighboring reports, transported to its frame. The old reports are complete on their input bodies at a join. On the new graph the first such update includes the cell's previous report, either from an old neighbor or, for a former singleton, from the cell's own new record. Unions on the fixed new graph are monotone with respect to report inclusion. Thus the reports only grow under any subset of updates, and every reported cell or bond is a true one. By fairness the reports stabilize. A fixed report contains all cell and bond records by connectedness. The cell recomputes its free labels from its report, so the settled labels and the permitted joins agree with the whole-body rule at every stock. The set of all finite reports is countable.

For the converse use independent replay of settled histories at every stock. Every two occurring bodies can be prepared side by side and must admit a join, and every body of a given mass can be prepared at its own stock. Thus the geometric catalogue is a merge-closed family, and the [whole-body reduction](rules.md#whole-body-reduction) builds it. ∎

## Simpler binding tables

A **matching table** gives each active label one possible partner, which can be itself. The local rule can then recover its own frozen label from the partner's facing-port label. The implementation proof works without reading its own bond record. Each new bond starts with matching labels, and recovery preserves those labels under any allowed order of updates.

**Corollary.** Matching tables suffice for every guaranteed target at a power-of-two stock. Take a recursive dissection from the [pairing theorem](fixed-target.md#pairing-and-powers-of-two). At each level choose one pair of edges witnessing the prescribed two-copy join. Label the two orbits of those edges in the current decorated body with a fresh pair of matching labels, or one self-matching label if the orbits coincide. Choosing other edges in the same orbits gives the same decorated outcome up to a proper motion. Disjoint label pairs at different levels allow only the prescribed pairings. The local implementation preserves those joins. ∎

A **self-binding table** allows labels to bind only themselves. This smaller class already loses targets at stock four. For example, a rule guaranteeing the four-square block must allow two identical reachable dominoes to join. A free label that can bind on a domino also allows its own copy to bind there, by a half-turn about the midpoint of that edge. This placement is disjoint. It makes a straight four at an end edge or an offset four on a side edge, neither the square. A matching table can instead pair complementary labels on the side edges. [Explicit four-square construction](repair.md#a-pulse-that-changes-the-next-repair)

**Open.** Whether matching tables can guarantee every whole-body target at every stock, and whether they realize every exact catalogue. Reproducing individual join histories is stronger. For example, at stock two a general table can admit numbered port pairs $0$–$1$, $1$–$2$, $2$–$3$ and refuse $0$–$3$. Under a matching table, the first two allowed contacts force ports 0 and 2 to have the same singleton label, making the fourth contact unavoidable. Both classes can still guarantee the unnumbered domino.

The theorem concerns assembly from identical loose cells. A comparison of writable memory also has to specify the allowed inputs and their duration. [A sharp example after a temporary input](repair.md#a-pulse-that-changes-the-next-repair)
