# Campaign manifests

Declarative descriptions of what to run: which cases, variants, sizes, builds
and run configurations, at how many samples, under which environment-quality
policy.

Two tiers (blueprint §15):

* **LEAN** — harness development, developer sanity, gross regression
  detection. Limited coverage and few samples. Every lean report must state
  prominently: **LEAN CAMPAIGN — NOT AUTHORITATIVE FOR CROSSOVER CLAIMS.**
* **FULL / AUTHORITATIVE** — all core workload families, practical size
  ladders, broad OpenMP scaling, CUDA SINGLE and DOUBLE where supported,
  sufficient independent samples, strict provenance and strict
  environment-quality gating.

A campaign manifest is a declaration, not a script: it names cells to execute
and leaves execution to the harness. The manifest format is defined by
**WP-09**; the first authoritative campaign is **WP-10**.

`campaign_id` from a manifest here is what appears in every record the
campaign produces.
