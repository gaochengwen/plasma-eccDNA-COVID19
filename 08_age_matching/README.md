# 08 — age-matched subset

**Figure S1** and **Table S4**.

| File | Role |
|---|---|
| `build_age_matched_pairs.py` | Construct the pair list: 1-year age caliper, no burden term. The pairing is written to `age_matched_pairs.json` and read from there by every consumer, so the figure and Table S4 can never diverge. |
| `rebuild_age_matched_pairs_v2_20260813.py` | Rebuild Table S4 on v2. Must run **after** Table S2a has been migrated to the v2 call set, because the matching inputs are read from it. |

Because matching uses age alone, the pair list is identical under v1 and v2; only the
paired burden values change.
