# Prefix-search results

The [certified runs](certified/) cover n = 3 through 16. They use
outward bounds in every geometric comparison. Each run reports zero
unmatched prefixes through stock 16.

The [prior numerical runs](prior-numerical/) give matching results
through n = 23. They use cyclotomic integer geometry and a
floating-point distance filter. The exact checker in this repository
can rerun any of those side counts. The longer runs can take tens of
minutes each.

The proof for all stock sizes also uses the
[continuation lemmas](../proofs/continuation.md).
