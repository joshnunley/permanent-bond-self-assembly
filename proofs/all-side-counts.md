# Classification beyond the fixed side counts

The [fixed-side-count checker](fixed-side-count.md) establishes the classification for $3\le n\le23$. A second proof takes all odd side counts at once. It uses the same continuation lemmas, with a search on turn words and short closures in place of a search at each $n$.

## Odd side counts

**Theorem.** For every odd $n\ge3$, the all-stock merge-closed families of regular $n$-gon bodies are precisely the admissible zigzag window families, with a single proper hand at each size, consistent up to global reflection.

**First five shapes.** Normalize headings to multiples of $\pi/n$. When $n$ is odd, the normals of one cell have a common parity and that parity flips across an edge contact. The contact graph is bipartite. A four-cycle would be a rhombus of unit side length. Its first and third directed edges would have opposite headings, which differ by $n$ in these units. But the two starting cells have the same parity. So cycles of length at most five are excluded.

At stock three the only connected tree is a path. At four the split $2+2$ excludes the three-leaf star. The possible nonpath trees at five are the four-leaf star and the tree with arms $1,1,2$. The star fails $2+3$. The other is a path of four with one leaf attached to its interior. The split $1+5$ at six adds a vertex to that tree. If it touches two old vertices it creates a cycle of length at most five, since the diameter of the old tree is three. Thus the new vertex is a leaf. The resulting types have arm lengths $(1,1,1,2)$, $(1,2,2)$, $(1,1,3)$, or two adjacent degree-three vertices. The first two and the last fail $3+3$ or $2+4$ respectively. The $(1,1,3)$ tree has a $2+4$ split only with a three-leaf star as the four-cell part, which is excluded. Hence the first five shapes are paths.

**From paths to the first nonpath.** The [chain-word checker](../code/chain_words_hands.py) describes a path by its sequence of signed turns. Proper congruence allows reading it backwards with all signs reversed. To grow a prefix, append a zero, a signed turn size already used, or a fresh turn size, at either end. The signs of fresh sizes are arbitrary names until their actual angles are assigned. For each extension try both hands and their union, and check all splits by cutting the word at the corresponding sizes. This enumeration covers all chain prefixes since a singleton extension of a path is at one end. At stock 14 the only surviving words are straight, constant-turn or alternating. The earlier shapes in such a prefix have that same turn pattern, since every smaller body can be a part of the stock-14 body by the merge condition.

A straight path through stock five forces a row, and an alternating path through stock fourteen forces a zigzag by the [continuation lemmas](continuation.md). A constant nonzero turn forces the same arc and eventually repeats a cell. If the first nonpath occurs through fourteen, it must be at a size from six to fourteen. The [late-case proof](late-cycles.md) checks these nine possibilities independently of $n$. None reaches stock fifteen. Thus an infinite odd-side family is a zigzag, since opposite normals required for a straight row are absent at odd $n$.

**Hands.** Odd windows are achiral. In any even window of at least four cells the first sign of its alternating turn word distinguishes its two hands, and reversing the traversal with sign reversal preserves that first sign. Two even pieces in a larger path must cut it at an even index. Their first signs agree. Pairing two different hands of the same size would therefore be impossible, and pairing any two even sizes forces their chosen hands to agree. The [zigzag construction](geometry.md#why-these-families-are-merge-closed) supplies the converse. ∎

As a finite consequence, for odd $n\ge7$ a sequence guaranteed through stock $2n+1$ already agrees with an all-stock zigzag family through that stock. At that length the early-tree argument applies and the late-case check excludes the first nonpath sizes 6–14. A constant turn cannot survive to $2n+1$, because every allowed heading is a multiple of $\pi/n$. The bound is sharp. Take a constant turn of $\pi/n$ to trace a regular $2n$-cell ring. Consecutive cells meet along whole edges. Nonconsecutive center distances are at least $2\cos(\pi/(2n))$ in units of the contact length, greater than the polygon's circumscribed diameter $\sec(\pi/n)$ for $n\ge7$. Every shorter stock can form the corresponding arc, with each split realized by consecutive pieces.

## Remaining even side counts

The general even-side extension has two early cases. The late-case computation, geometric continuation and realizability of the listed families are already available.

1. **First nonpath at stock three or four.** Certify that a surviving prefix agrees with a two-row ladder far enough to apply continuation.
2. **First nonpath at stock five.** Certify that every candidate fails to continue to the required cutoff.

Those early angle searches have numerical evidence. Completing this extension requires certifying the search over angle intervals and exceptional parameter values, or a replacement proof.
