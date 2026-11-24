# Determinate assembly with permanent bonds

Start with any number of identical regular polygons. Two rigid bodies
can join along a whole edge when the two labels are compatible. Labels
can depend on the structure visible from the joining edge. Can one
rule, used at every stock size, force every maximal assembly history to
finish with the same shape for that stock size?

For regular polygons with **3 through 16 sides**, the answer is the
same in all three label models. A final shape is an arrangement of cells
up to rotation, translation and reflection.

| Information available to each label | Any fixed finite radius | Whole body |
| --- | --- | --- |
| Whole-body edge numbering (A) | Straight row for even side counts. No determinate rule for odd side counts. | Rows, zigzag windows and two-row ladders |
| Only the marked edge numbered (C) | Same | Same |
| Blank cells (B) | Same | Same |

The whole-body list consists of straight rows for even side counts,
contiguous windows of alternating zigzags with bond separation
n/6 ≤ k < n/2, and two-row ladders for even side counts with a rung
angle between 60° and 90° inclusive.
The two ladder variants differ in which row gets the extra cell when
the stock is odd.

![Nine decagons arranged in a row, zigzag, and the two ladder variants](figures/families.png)

*At stock nine, the odd ladder cell is highlighted.*

The **finite-radius obstruction** holds for every side count. No
nonstraight zigzag window family or two-row ladder family can be built
by a determinate rule at any fixed finite radius, even with the full
numbering visible in that radius. A blank radius-zero rule builds the
straight row whenever the side count is even. The **whole-body
classification** uses a finite exact computation for 3 through 16 sides
and a proof for all stock sizes. Extending the classification to every
side count still requires certifying two angle-parameter cases. Earlier
searches also matched the list for 17 through 23 sides. Their outward-rounded reruns
are tracked in [the results](results/).

## Proof and reproduction

- [The three label models and the reduction from rules to geometry](proofs/rules.md)
- [Definitions of rows, zigzag windows and ladders, with their constructions](proofs/geometry.md)
- [Exhaustive prefix check for each side count](proofs/fixed-side-count.md)
- [Continuation past the finite cutoff and consistency of hands](proofs/continuation.md)
- [Finite-radius obstruction and the radius-zero row rule](proofs/finite-radius.md)
- [Status of the extension to all side counts](proofs/all-side-counts.md)

From this directory, with Python 3.11 or later, install the dependencies
and run an example check.

```sh
python3 -m pip install -r requirements.txt
python3 -B code/exact_n.py 8 16 --require-match
```

To rerun the classification, run the same command for each side count.

```sh
for n in $(seq 3 16); do
  python3 -B code/exact_n.py "$n" 16 --require-match || exit 1
done
```

Each check should report `unmatched prefixes: 0`. A listed constant-turn
arc can occur in a finite prefix. The continuation proof excludes it
from an infinite family. See [the recorded runs](results/).

## Related models

The polygonal free-body tile assembly model also permits rotatable
polygonal tiles. The question here fixes regular polygons and asks
which shape families one unseeded rule can force at every finite
stock, across every assembly history.
[Demaine et al., *One Tile to Rule Them All*](https://arxiv.org/abs/1212.4756)

Pairwise joining of independently formed bodies is shared with
hierarchical tile assembly. In this model, permanence and exact
finite stock imply the merge condition on every split of a stock.
[Doty, *Producibility in Hierarchical Self-Assembly*](https://arxiv.org/abs/1304.7804)
