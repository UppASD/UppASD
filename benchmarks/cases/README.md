# Benchmark cases

Immutable case definitions: one directory per case family, each holding a case
manifest and the input templates it owns.

Planned families (blueprint §3):

| Case | Workload |
| --- | --- |
| `B01_bccFe` | bcc Fe, medium-range interactions — central reference case |
| `B02_skyrmion2D` | 2D skyrmion, short-range J+D |
| `B03_skyrmion3D` | 3D skyrmion/chiral magnet, short-range J+D |
| `B04_dhcpNd` | dhcp Nd, very long-range interactions |
| `B05_dipoleFFT` | dipole/FFT workload |

Rules that already apply:

* **Templates are immutable inputs.** A benchmark run never rewrites anything
  here; it copies what it needs into its own generated work directory.
* **No generated output in this directory.** Results and work directories live
  under gitignored paths — see the result-data policy in [`../README.md`](../README.md).
* **No physics simplification for convenience.** Interaction shells are not
  truncated, cutoffs are not reduced, weak interactions are not discarded and
  the basis is not simplified in order to make benchmarking easier.

The case manifest format — variants, legal size ladders, workload metadata and
input override allow-lists — is defined by **WP-02**; individual families are
admitted by **WP-08A**…**WP-08E**.
