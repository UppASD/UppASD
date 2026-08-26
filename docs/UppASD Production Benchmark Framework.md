# UppASD Production Benchmark Framework

## Master Blueprint and Governance Document

### Repository scope

UppASD production benchmarking framework.

Primary execution modes considered in this project:

* CPU/OpenMP
* CUDA
* HIP only where an implemented backend exists and suitable hardware is available

Adaptive coarse graining may later be benchmarked as an additional UppASD workflow using the same framework, but it is **not** the organizing principle of this benchmark project.

---

# 0. Mandatory MPI exclusion

> **Disregard all MPI-related statements in this document, in the associated work-package prompts, and in any derived benchmark plans or actions. MPI is not implemented in UppASD for the scope considered here and is neither to be pursued nor benchmarked as part of this project.**

This instruction overrides any residual MPI wording that may remain in earlier versions of the benchmark blueprint or associated WP prompts.

Specifically:

* do not add MPI as a benchmark dimension;
* do not add `mpi_ranks` or equivalent fields to schemas;
* do not add MPI launcher support;
* do not use `mpirun`, `mpiexec`, or equivalent launch mechanisms;
* do not reserve infrastructure specifically for future MPI benchmarking;
* do not create MPI scaling campaigns;
* do not modify UppASD to add MPI support;
* do not include MPI as deferred work for this project.

The CPU performance model considered here is **shared-memory CPU/OpenMP only**.

---

# 1. Purpose

Build a general, reproducible production-performance framework for UppASD that answers practical questions about:

* CPU/OpenMP scaling;
* GPU scaling;
* CPU ↔ GPU crossover;
* SINGLE versus DOUBLE precision behaviour;
* problem-size scaling;
* interaction-range and connectivity dependence;
* neighbour-list versus FFT/dipolar workloads;
* setup cost versus steady-state simulation cost;
* measurement overhead;
* machine-to-machine performance differences.

The framework should ultimately provide a trustworthy basis for answering questions such as:

1. How many OpenMP threads should UppASD use for a given workload?
2. At what problem size does CUDA become faster than the best CPU/OpenMP configuration?
3. How does that crossover depend on interaction range?
4. How different are SINGLE and DOUBLE precision GPU performance?
5. Is atom count or interaction count the more useful predictor of performance?
6. How does dipole/FFT performance differ from ordinary short-range Hamiltonian evaluation?
7. Which production components actually deserve further optimization?

---

# 2. Core benchmarking principle

Headline performance numbers must come from the **real UppASD production executable and complete production workflows**.

Kernel benchmarks, profiler measurements, isolated microbenchmarks, and phase timings are useful secondary diagnostic evidence, but they must never replace the complete production timing in claims such as:

* CPU speedup;
* GPU speedup;
* CPU/GPU crossover;
* precision crossover;
* workflow crossover.

The numerator and denominator of every reported speedup must represent equivalent physical and numerical work.

For ordinary ASD, the complete production timestep includes the appropriate production sequence, conceptually:

[
H(\mathbf m_n)
\rightarrow
\text{Depondt predictor}
\rightarrow
H(\mathbf m^*)
\rightarrow
\text{Depondt corrector},
]

including finite-temperature Langevin/thermal-field work whenever the selected production configuration uses it.

Do not derive a production speedup by comparing a complete timestep on one backend against a subset of the timestep on another.

---

# 3. Initial benchmark workload families

The initial benchmark suite consists of five workload families chosen to span distinct computational regimes.

## B01 — bcc Fe: medium-range interactions

Purpose:

* canonical three-dimensional ASD workload;
* realistic medium-range exchange interaction cloud;
* primary OpenMP scaling reference;
* primary CPU/GPU crossover reference;
* finite-temperature Depondt/Langevin benchmark.

Suggested variants:

* (T=0);
* finite-temperature production dynamics;
* optional controlled measurement cadence.

This should normally serve as the central reference case.

---

## B02 — 2D skyrmion: short-range (J+D)

Purpose:

* short neighbour lists;
* exchange + DMI;
* quasi-two-dimensional geometry;
* representative chiral/topological-magnet workload;
* comparison against three-dimensional short-range dynamics.

The benchmark significance is the computational workload, not the exact skyrmion texture.

---

## B03 — 3D skyrmion/chiral magnet: short-range (J+D)

Purpose:

* genuinely three-dimensional short-range (J+D) dynamics;
* contrast with B02 while retaining a similar Hamiltonian class;
* representative three-dimensional chiral-spin calculation;
* later useful as a realistic adaptive-CG workload.

A separate Hopfion benchmark is intentionally omitted from the initial benchmark suite because its ordinary ASD computational workload is sufficiently represented by this three-dimensional short-range (J+D) family.

---

## B04 — dhcp Nd: very long-range interactions

Purpose:

* deliberately large interaction cloud;
* neighbour traversal stress;
* memory-bandwidth/cache stress;
* OpenMP saturation study;
* GPU work-amortization study;
* test of whether total interaction count predicts CPU/GPU crossover better than atom count.

The physically intended long-range interaction cloud must be retained.

Do **not**:

* truncate neighbour shells;
* reduce exchange cutoff;
* discard weak interactions;
* simplify the basis;

merely to make benchmarking easier.

This workload is important regardless of its suitability for adaptive coarse graining.

---

## B05 — dipole/FFT workload

Purpose:

* represent an algorithmically different production workload;
* characterize FFT-based dipole/magnetostatic performance;
* quantify FFT setup and steady-state costs separately;
* determine CPU/GPU crossover for a global rather than purely neighbour-list calculation.

Prefer a physically meaningful but reasonably simple model in which the dipolar workload is visible rather than hidden beneath unrelated expensive Hamiltonian terms.

Record FFT-grid information explicitly.

---

# 4. Benchmark dimensions

The framework must support independent variation of:

* benchmark family;
* physical variant;
* system size;
* backend;
* precision;
* OpenMP thread count;
* measurement profile;
* simulation length;
* machine/hardware configuration.

Initial production backends:

* CPU/OpenMP
* CUDA

HIP should be represented only to the extent that UppASD already provides an applicable implementation. No HIP performance claim is to be made without actual suitable hardware execution.

Initial precision variants:

* SINGLE
* DOUBLE

`MIXED` remains unsupported while the production build deliberately does not implement it.

Do not implement mixed precision as part of this benchmark project.

---

# 5. Workload size is multidimensional

Atom count alone is not a sufficient description of UppASD computational work.

For every neighbour-list benchmark, record at minimum:

[
N_{\rm atom},
]

[
N_{\rm directed\ interactions},
]

[
\langle N_{\rm neigh}\rangle,
]

and

[
N_{\rm neigh}^{\rm max}.
]

Where inexpensive, also record the median neighbour count.

The principal neighbour-workload analyses should therefore include both:

[
T(N_{\rm atom})
]

and

[
T(N_{\rm directed\ interactions}).
]

This is especially important when comparing short-range skyrmion systems with long-range dhcp Nd.

## FFT/dipole workloads

For dipole/FFT cases additionally record:

* FFT grid dimensions;
* total FFT grid points;
* relevant padding or grid-expansion metadata where available;
* dipole implementation/backend.

Analyze these workloads against both atom count and FFT-grid size.

Do not apply neighbour-interaction throughput metrics blindly to FFT workloads.

---

# 6. Timing semantics

Every benchmark should distinguish three performance regimes.

## 6.1 Full process time

[
T_{\rm process}
]

Includes:

* executable startup;
* input parsing;
* allocations;
* Hamiltonian setup;
* neighbour setup;
* GPU initialization;
* FFT planning where applicable;
* simulation;
* requested measurements;
* finalization.

This quantity determines the economics of short jobs.

---

## 6.2 Setup/fixed-cost estimate

Estimate the approximately fixed component separately.

The preferred initial model is:

[
T(n)
====

T_{\rm fixed}
+
n,t_{\rm step}.
]

The fitted intercept should initially be called something such as:

`setup_or_fixed_intercept`

rather than being interpreted automatically as pure setup cost.

---

## 6.3 Steady-state production timestep

[
t_{\rm step}
]

This is the primary compute-performance quantity.

The initial preferred measurement method is executable-level multi-(N_{\rm step}) timing using several complete production runs.

For authoritative campaigns, use at least three sufficiently separated step counts.

A lean developer tier may use a calibrated two-point subtraction.

Do not construct steady-state performance by subtracting unrelated profiler phase timings.

---

# 7. Measurement profiles

At minimum provide two explicitly different profiles.

## DYNAMICS_ONLY

Suppress optional measurement/output work as far as normal UppASD input semantics permit.

This is the primary computational benchmark.

It should still execute all work inherently required by the chosen dynamics and Hamiltonian.

---

## PRODUCTION_LIGHT

Use a fixed, documented, representative measurement cadence.

Purpose:

* quantify realistic measurement overhead;
* understand when measurement/output becomes important;
* estimate real simulation throughput.

Do not mix `DYNAMICS_ONLY` and `PRODUCTION_LIGHT` results in the same crossover curve without explicit labelling.

---

# 8. CPU/OpenMP baseline methodology

CPU performance must not be represented by one arbitrary OpenMP thread count.

For every workload and problem size measure:

[
T(N,p),
]

where (p) is the OpenMP thread count.

Calculate:

[
S_{\rm OMP}(p)
==============

\frac{T(N,1)}{T(N,p)}
]

and

[
E_{\rm OMP}(p)
==============

\frac{S_{\rm OMP}(p)}{p}.
]

Determine experimentally:

[
p_{\rm best}(N,\mathrm{case}),
]

the fastest valid measured CPU/OpenMP configuration.

Maintain two canonical CPU references.

## CPU-1T

One OpenMP thread.

Useful for:

* historical comparisons;
* algorithmic comparisons;
* assessing parallel efficiency.

## CPU-BEST

The fastest measured supported OpenMP configuration for that case and size.

This is the principal CPU baseline for GPU crossover.

A GPU speedup against one CPU thread may be reported, but must not replace GPU/CPU-BEST as the production comparison.

---

# 9. OpenMP execution control

CPU benchmarking must use explicit and reproducible OpenMP configuration.

Record at minimum:

* `OMP_NUM_THREADS`;
* `OMP_PLACES`;
* `OMP_PROC_BIND`;
* `OMP_DYNAMIC`;
* process affinity;
* CPU topology;
* NUMA topology.

Prioritize physical cores in the standard campaign.

SMT/hyperthreading may be investigated separately but must not be mixed silently with physical-core scaling.

Multi-socket results must identify whether execution spans:

* one socket/NUMA region;
* several NUMA regions;
* the full node.

---

# 10. GPU crossover

For every relevant case/configuration calculate:

[
S_{\rm GPU/BESTCPU}
===================

\frac{T_{\rm CPU-BEST}}
{T_{\rm GPU}}.
]

Also retain:

[
S_{\rm GPU/1T},
]

but treat it as a secondary metric.

Determine where measured/interpolated performance crosses:

* (1.0\times);
* (1.25\times);
* (2.0\times);
* (5.0\times).

Interpolation is allowed only between valid measured neighbouring points.

Never extrapolate outside the measured range and present the result as a measured crossover.

Use classifications such as:

* `below_tested_range`;
* `within_tested_range`;
* `above_tested_range`.

---

# 11. Precision semantics

The benchmark framework must distinguish:

* requested precision;
* effective CPU precision;
* effective GPU precision;
* relevant measurement/reduction precision.

Do not infer numerical precision solely from CMake option names.

Before precision comparisons are accepted, audit what the production configuration actually changes.

Classify comparisons as either:

## PRECISION_MATCHED

The relevant CPU and GPU numerical paths use genuinely corresponding precision.

## PRODUCTION_CONFIGURATION

The comparison represents real supported production modes but not identical precision.

Example:

* CPU double;
* CUDA single.

Such a comparison can be operationally useful, but it must not be labelled precision-matched.

---

# 12. Throughput metrics

For all workloads report:

[
R_{\rm spin}
============

\frac{N_{\rm atom}}{t_{\rm step}}
]

as spin-steps/s or an equivalently clear name.

For neighbour-list workloads additionally report:

[
R_{\rm interaction}
===================

\frac{N_{\rm directed\ interactions}}
{t_{\rm step}}.
]

Use a name such as:

`directed_interaction_visits_per_second`

unless source inspection supports a more exact operation interpretation.

For FFT/dipole workloads report a suitable grid-normalized metric such as:

[
R_{\rm FFTgrid}
===============

\frac{N_{\rm FFT\ grid\ points}}
{t_{\rm step}}.
]

Do not force every workload into one normalization.

---

# 13. Performance provenance

Every authoritative result must identify, directly or through a referenced immutable build/campaign record:

## Source

* git commit;
* tracked clean/dirty state;
* changed tracked files if dirty;
* binary checksum.

## Build

* compiler;
* compiler version;
* relevant compile flags;
* CMake configuration;
* backend;
* requested precision;
* effective precision.

## CPU

* model;
* sockets;
* physical cores;
* SMT state;
* NUMA topology;
* affinity;
* OpenMP environment.

## GPU

Where applicable:

* model;
* device ID/UUID;
* driver;
* CUDA/HIP runtime/toolkit;
* memory;
* temperature;
* clocks;
* utilization;
* throttling status;
* evidence of competing processes where obtainable.

## Workload

* case;
* variant;
* size;
* atom count;
* interaction count or FFT-grid size;
* temperature;
* timestep;
* measurement profile.

## Statistics

* raw individual timing samples;
* sample count;
* median;
* MAD;
* minimum;
* maximum.

Do not report only the fastest sample.

---

# 14. Environment quality

Separate:

`numerical_valid`

from:

`environment_valid`.

Examples:

A simulation that completes correctly while another process occupies the GPU is:

* numerically valid;
* environmentally unsuitable for authoritative timing.

Possible environment-quality flags include:

* `dirty_tree`;
* `gpu_busy`;
* `gpu_thermal_throttle`;
* `gpu_clock_unstable`;
* `cpu_affinity_unknown`;
* `cpu_frequency_unstable`;
* `background_load_high`;
* `metadata_incomplete`.

Authoritative campaigns should support a strict environment mode that refuses unsuitable samples.

Developer runs should remain usable even when some quality flags are present.

---

# 15. Campaign tiers

## LEAN

Purpose:

* harness development;
* developer sanity;
* gross regression detection;
* rapid local comparisons.

Use:

* a limited case subset;
* few sizes;
* limited OpenMP thread counts;
* small sample count;
* selected GPU configurations.

Every lean report must state prominently:

> **LEAN CAMPAIGN — NOT AUTHORITATIVE FOR CROSSOVER CLAIMS**

---

## FULL / AUTHORITATIVE

Purpose:

* production OpenMP characterization;
* CPU/GPU crossover;
* precision effects;
* interaction-range analysis;
* dipole/FFT analysis;
* release/performance documentation.

Requirements should include:

* all five core workload families;
* practical size ladders;
* broad OpenMP scaling;
* CUDA SINGLE and DOUBLE where supported;
* sufficient independent samples;
* strict provenance;
* strict environment-quality gating.

---

# 16. Continuous integration policy

Do not impose fixed performance thresholds on generic/shared CI hardware.

Normal CI may verify:

* schema parsing;
* case-manifest validation;
* benchmark-runner logic;
* result aggregation;
* report generation;
* tiny infrastructure-only benchmark fixtures.

Authoritative performance regression testing requires controlled hardware and is outside ordinary shared CI timing.

---

# 17. Repository structure

Preferred structure:

```text
benchmarks/
    README.md
    schema/
    cases/
    campaigns/
    harness/
    analysis/
    tests/
```

Generated runtime data should normally live in:

* gitignored result directories; or
* dedicated external campaign result directories.

Production benchmark templates are immutable inputs.

Every executable run uses its own generated work directory.

Benchmark execution must never rewrite tracked production templates.

---

# 18. Work-package structure

The project is divided into the following work packages.

## Infrastructure

**WP-01 — Benchmark contract and result schema**

Establish:

* repository skeleton;
* terminology;
* machine-readable schema;
* validity concepts.

---

**WP-02 — Declarative case/template infrastructure**

Establish:

* immutable templates;
* case manifests;
* variants;
* legal size scaling;
* workload metadata;
* input override allow-lists.

---

**WP-03 — Production executable timing**

Establish:

* isolated runs;
* complete executable timing;
* multi-(N_{\rm step}) steady-state estimation;
* dynamics-only and production-light profiles.

---

**WP-04 — Provenance and environment quality**

Establish:

* source/build provenance;
* CPU/GPU metadata;
* contention/throttling detection;
* sample-quality classification.

---

**WP-05 — OpenMP scaling**

Establish:

* thread sweeps;
* affinity/binding;
* CPU-1T;
* CPU-BEST;
* OpenMP efficiency.

---

**WP-06 — GPU and precision campaigns**

Establish:

* CUDA SINGLE/DOUBLE;
* precision audit;
* matched versus production-configuration comparisons;
* GPU memory-limit handling.

---

**WP-07 — Scaling/crossover analysis**

Establish:

* crossover detection;
* atom-normalized throughput;
* interaction-normalized throughput;
* FFT-grid analysis;
* automated plots/reports.

---

## Production workload admission

Once WP-02 is stable, the benchmark families can be admitted independently.

**WP-08A — bcc Fe medium-range**

**WP-08B — 2D skyrmion short-range (J+D)**

**WP-08C — 3D skyrmion/chiral short-range (J+D)**

**WP-08D — dhcp Nd long-range**

**WP-08E — dipole/FFT**

---

## Campaign execution

**WP-09 — Lean and full campaign definitions**

Create standardized campaign manifests.

---

**WP-10 — First authoritative production characterization**

Run the complete benchmark campaign and establish:

* OpenMP scaling;
* CPU-BEST;
* CUDA crossover;
* precision effects;
* interaction-range dependence;
* FFT/dipole behaviour;
* setup versus steady-state economics;
* evidence-based optimization priorities.

---

# 19. Dependency graph

```text
MASTER BLUEPRINT
       |
       v
     WP-01
       |
       v
     WP-02
       |
       +---------------------> WP-08A  bcc Fe
       +---------------------> WP-08B  2D skyrmion
       +---------------------> WP-08C  3D skyrmion
       +---------------------> WP-08D  dhcp Nd
       +---------------------> WP-08E  dipole/FFT
       |
       v
     WP-03
       |
       v
     WP-04
       |
       v
     WP-05
       |
       v
     WP-06
       |
       v
     WP-07
       |
       +---- requires production families before authoritative campaign
       |
       v
     WP-09
       |
       v
     WP-10
```

WP-08A through WP-08E may run independently once WP-02 has frozen the case/template contract.

WP-03 through WP-07 should normally remain sequential because each establishes assumptions consumed by the next.

WP-10 must not start until:

* production timing semantics are accepted;
* provenance/environment gates are accepted;
* all core workload families required by the campaign are admitted.

---

# 20. Explicit initial-project exclusions

The following are **not** part of this project:

* MPI support;
* MPI implementation;
* MPI benchmarking;
* mixed-precision implementation;
* adaptive-CG-specific optimization;
* production kernel optimization;
* speculative HIP performance claims without hardware;
* artificial simplification of benchmark Hamiltonians;
* hard performance gating on uncontrolled shared CI machines.

If later work adds new execution paradigms to UppASD, benchmark support can be considered at that time rather than pre-designed here.

---

# 21. Adaptive coarse graining

Adaptive coarse graining should eventually enter the framework as another production workflow, for example conceptually:

```text
workflow = ordinary_asd
workflow = adaptive_cg
```

Both must obey exactly the same headline timing contract.

Additional adaptive dimensions may include:

* block size;
* selector configuration;
* actual effective fine fraction;
* interface fraction;
* transition cadence.

Forced fine fractions may be useful diagnostic controls but must be clearly distinguished from naturally selected production adaptive states.

The general benchmark infrastructure must not be modified to make adaptive CG appear favourable.

---

# 22. Performance optimization policy

Do not perform production optimization during construction of the benchmark framework unless required to make benchmarking function correctly.

The first authoritative campaign should benchmark the current production implementation.

If a bottleneck is discovered:

1. record it;
2. reproduce it;
3. freeze the benchmark evidence;
4. then create a separate optimization task;
5. use the same benchmark framework to measure the improvement.

This avoids optimizing toward unstable or incomplete microbenchmarks.

---

# 23. Success criteria

The framework is considered successful when it can reproducibly answer:

1. How does CPU/OpenMP scaling depend on workload and system size?
2. What OpenMP configuration gives CPU-BEST for each case?
3. At what size does CUDA outperform CPU-BEST?
4. At what size does CUDA reach useful speedups such as (2\times) and (5\times)?
5. How does SINGLE versus DOUBLE change the crossover?
6. How does short-, medium-, and long-range interaction structure affect performance?
7. Is total directed-interaction count a better workload predictor than atom count?
8. How does dipole/FFT scaling differ from local neighbour-list Hamiltonians?
9. When does initialization/setup cost matter?
10. What production components should be optimized next?

The final purpose of the framework is not merely to create benchmark plots.

It is to provide reliable evidence for future UppASD engineering decisions.
