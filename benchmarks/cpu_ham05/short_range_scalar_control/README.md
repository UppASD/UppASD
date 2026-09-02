# CPU-HAM-05 short-range scalar-J control

This campaign-only fixture is a one-atom simple-cubic lattice with six
nearest-neighbour scalar-J bonds, fully periodic in all three dimensions.
It exists to cover the low-coordination regime in the CPU-HAM-05 crossover
matrix without changing the five admitted general benchmark families.

The fixture is deterministic (`temp=0`, `ip_mode=N`) and has no DMI,
anisotropy, dipole, or other pair term. The campaign runner changes only
`ncell`, `nstep`, measurement cadence, and explicit CPU backend-dispatch
keywords in generated copies; these files are never run in place.
