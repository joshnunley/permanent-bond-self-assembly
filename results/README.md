# Reproducing the results

Use Python 3.11 or later from the repository root. The geometric verifier uses NumPy, SymPy and mpmath.

```sh
python3 -m pip install -r requirements.txt
python3 -B code/exact_n.py 8 16 --require-match
```

## Fixed side count

The [certified runs](certified/) use outward bounds for geometric comparisons. Each reports zero unmatched prefixes through stock 16. [Proof and search design](../proofs/fixed-side-count.md)

```sh
for n in $(seq 3 23); do
  python3 -B code/exact_n.py "$n" 16 --require-match || exit 1
done
```

Runs at the larger side counts can take tens of minutes. The [prior numerical runs](prior-numerical/) are retained as historical records.

## Arbitrary odd side counts

The [odd-side proof](../proofs/all-side-counts.md#odd-side-counts) uses three finite checks whose input is a stock cutoff or a candidate first nonpath size. It does not loop over side counts.

```sh
python3 -B code/chain_words_hands.py 14
python3 -B code/first_nonchain_graphs.py 14
python3 -B code/late_cycle_extensions.py 6 7 8 9 10 11 12 13 14
```

The first check has only constant, alternating and zero turn words at stock 14. The second checks the possible contact graphs of a first nonpath. The third excludes continuation to 15 when the first nonpath occurs at any size from 6 to 14. [Chain-word record](one_turn_words.txt) · [Contact-graph record](first_nonchain_graphs_14.txt) · [Late-case record](late_cycle_extensions.txt)

## Fixed-stock atlas

The atlas has its own [reproduction instructions](../atlas/README.md#reproduction), using only the Python standard library.
