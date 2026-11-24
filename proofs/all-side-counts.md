# Extension to every side count

The target is the [classification theorem](fixed-side-count.md) with
every n ≥ 3 in place of the certified interval. The families on its
list are realizable for every admissible side count. The finite-radius
exclusions and the continuation arguments also hold for every side count.

For the remaining necessity argument, split by the first stock size
whose contact graph is neither a path nor a singleton. The proof draft uses a chain-word
reduction to bound that first size by 14. Its short-closure enumeration
and extension checks cover first sizes 6 through 14. In particular, correcting a sign in the relation
between a short closure and its angles changes the cases for first sizes
6 and 7. The corrected checks both still exclude continuation to stock 15.

The remaining certification work consists of two early cases in the
angle-parameter analysis. Both have positive numerical evidence.

1. **First non-chain size 3 or 4.** Certify that a surviving prefix
   agrees with a two-row ladder far enough to apply continuation.
2. **First non-chain size 5.** Certify that every candidate fails to
   continue to the required cutoff.

These two angle searches currently use numerical screens for angle boxes
and exceptional roots. The all-side-count classification requires a
completeness proof for those screens or an exact replacement. The
[fixed-side-count check](fixed-side-count.md) uses a different search
that runs directly in a fixed cyclotomic field, so it avoids those
parameter screens.
