# Campaign manifests

Declarative descriptions of what to run: which cases, variants, sizes,
backends/precisions and run configurations, at how many samples, under which
environment-quality policy. Schema:
[`../schema/campaign_manifest.v1.schema.json`](../schema/campaign_manifest.v1.schema.json).
Loaded and cell-resolved by
[`../harness/campaigns.py`](../harness/campaigns.py). Neither this directory
nor that module executes anything -- a campaign manifest names cells; a
campaign driver (WP-10) hands them to `harness.runner`/`harness.gpu_campaign`.

Two tiers (blueprint §15, TERMINOLOGY.md CAMPAIGN):

* **LEAN** — harness development, developer sanity, gross regression
  detection. Limited coverage and few samples. Every lean report must state
  prominently: **LEAN CAMPAIGN — NOT AUTHORITATIVE FOR CROSSOVER CLAIMS.**
  Enforced by the schema (a `LEAN` manifest requires `report_banner` set to
  exactly that text) and rendered as the first line of any report built with
  `analysis.markdown_report.generate_markdown_report(..., report_banner=...)`.
* **FULL / AUTHORITATIVE** — all core workload families, practical size
  ladders, broad OpenMP scaling, CUDA SINGLE and DOUBLE where supported,
  at least 5 samples per cell, strict provenance and strict
  environment-quality gating. A `FULL` manifest may never carry
  `report_banner`, may never resolve to an `infrastructure_test_only` case
  (`harness.cases.require_not_infrastructure_only`), and must set
  `environment_quality_mode: STRICT`.

`campaign_id` from a manifest here is what appears in every record the
campaign produces.

## The four standard campaigns

| File | Tier | Backend scope | Purpose |
| --- | --- | --- | --- |
| [`lean.yaml`](lean.yaml) | LEAN | CPU + CUDA | Harness sanity / gross regression detection. bcc Fe + one short-range J+D case, 3 sizes, a 3-point OpenMP sample, CUDA SINGLE+DOUBLE. |
| [`full_cpu.yaml`](full_cpu.yaml) | FULL | CPU only | All 5 families, full size ladders, broad physical-core OpenMP sweep -- establishes CPU-1T and CPU-BEST. |
| [`full_cuda.yaml`](full_cuda.yaml) | FULL | CUDA only | All 5 families where GPU memory permits, SINGLE + DOUBLE, strict GPU environment gate. |
| [`full_crossover.yaml`](full_crossover.yaml) | FULL | CPU + CUDA | **The manifest WP-10 actually runs.** CPU-BEST and CUDA together under one `campaign_id`/build/source revision (`single_campaign_execution: true`), so every crossover computed from its aggregates is a real same-run comparison. |

`full_cpu.yaml` and `full_cuda.yaml` are useful independently (e.g. a
CPU-only OpenMP-scaling study, or a CUDA-only precision study), but their
records must never be paired against each other for a crossover claim
(blueprint §10: "Do not rerun an arbitrary different CPU configuration
merely because the GPU campaign was started separately") --
`analysis.crossover.compute_gpu_speedup` already refuses to compare two
aggregates that do not share one `campaign_id`, which structurally forces
crossover claims onto `full_crossover.yaml`'s own records.

[`ci_infra.yaml`](ci_infra.yaml) is a fifth, deliberately tiny manifest:
`INFRA_test_only` only, LEAN tier, CPU only. It exists solely so ordinary
shared CI can exercise manifest parsing, cell resolution and the runner
(blueprint §16) without a wall-time pass/fail threshold and without ever
being mistaken for benchmark evidence.

## Machine hardware policy is kept separate from case selection

A campaign manifest's `cpu` block never names a concrete thread-count list.
It names a `thread_policy` -- `LEAN_SAMPLE` (1 / moderate / full
physical-core count) or `FULL_PHYSICAL_SWEEP` (every count from 1 to the
host's physical-core count) -- and `harness.campaigns.resolve_cpu_sweep`
turns that into a real `harness.omp_sweep.OmpSweep` only once a real
`harness.omp_topology.HostTopology` is available. The same manifest is
therefore portable across any machine; no host's core count is ever checked
into a tracked campaign file.
