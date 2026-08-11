"""Source-of-truth fixture names shared by the adaptive-CG harnesses."""

FEATURE_OFF_CASE = "feature_off"
STATIC_ALL_FINE_CASE = "static_all_fine"
STATIC_ALL_COARSE_CASE = "static_all_coarse"
STATIC_MIXED_CASE = "static_mixed"
ADAPTIVE_MIXED_CASE = "adaptive_mixed"
POLARIZATION_GATE_CASE = "polarization_gate_cpu"
POLARIZATION_GATE_GPU_CASE = "polarization_gate_gpu"

CPU_REJECTION_CASES = (
    "invalid_partial_block",
    "invalid_mask",
    "unsupported_temperature",
    "unsupported_initial_phase_x",
    "missing_alat",
)
CPU_ANISOTROPY_CASES = (
    "dmi_anisotropy_mixed",
    "anisotropy_uniform_fine",
    "anisotropy_uniform_coarse",
)
# RCG-03 ANI-NONUNIFORM-REJECT: a do_cluster embedding placed exactly on an
# existing host lattice site (clus_expand=0, so the pre-existing
# Natom/=NA*N1*N2*N3 geometry rejection never triggers) that gives one atom
# a divergent anisotropy from every other atom sharing its basis index.
CLUSTER_ANISOTROPY_REJECTION_CASE = "ani_nonuniform_reject_cluster"
CPU_PARITY_CASES = {"static": "parity_static_cpu", "adaptive": "parity_adaptive_cpu"}
CPU_INITIAL_PHASE_CASES = {
    "M": ("initial_phase_mc", "Performing MC initial phase:"),
    "H": ("initial_phase_heat_bath", "Performing MC initial phase:"),
    "Q": ("initial_phase_q", "Line search minimization done"),
    "Y": ("initial_phase_y", "Line search minimization done"),
    "Z": ("initial_phase_z", "Line search minimization done"),
    "G": ("initial_phase_g", "Performing initial phase: energy minimization"),
}
GPU_EXECUTABLE_CASES = ("gpu_static_mixed", "gpu_adaptive_mixed")
GPU_INITIAL_PHASE_CASES = (
    "initial_phase_sd_gpu", "initial_phase_mc_gpu", "initial_phase_q_gpu"
)
GPU_PARITY_CASES = {"static": "parity_static_gpu", "adaptive": "parity_adaptive_gpu"}
GPU_DMI_CASE = "dmi_anisotropy_mixed_gpu"
GPU_FFT_CASE = "gpu_fft_static_mixed"
CPU_INITIAL_PHASE_SD_CASE = "initial_phase_sd_cpu"
SPIN_SPIRAL_CASE = "initmag_spin_spiral"
RESTART_INITMAG_CASE = "initmag_restart_atomistic"
EXAMPLE_CASES = ("static_mixed", "adaptive", "initial_phase_texture")

# RCG-04D E2E-MOVING-OFF-FINE: genuinely nonstationary (nonzero-torque)
# conical-spiral feature-off/all-fine parity pair, distinct from the
# uniform/stationary FEATURE_OFF_CASE/STATIC_ALL_FINE_CASE pair above (which
# remains a zero-torque smoke test only). See run_moving_off_fine.py and
# each fixture's README.md.
MOVING_FEATURE_OFF_CASE = "moving_feature_off"
MOVING_ALL_FINE_CASE = "moving_all_fine"

# RCG-04E E2E-MOVING-ALL-COARSE: wide-geometry (ncell 24 2 2) re-run of the
# RCG-04D feature-off/all-fine parity pair, used to (re-)establish the
# atomistic oracle at this slice's geometry, plus four all-coarse block-size
# resolutions compared against it. See run_moving_all_coarse.py and each
# fixture's README.md.
MOVING_FEATURE_OFF_WIDE_CASE = "moving_feature_off_wide"
MOVING_ALL_FINE_WIDE_CASE = "moving_all_fine_wide"
MOVING_ALL_COARSE_CASES = (
    "moving_all_coarse_bs1", "moving_all_coarse_bs2",
    "moving_all_coarse_bs4", "moving_all_coarse_bs8",
)

# RCG-04F E2E-MOVING-STATIC: static fine/interface(buffer)/coarse mixed-mask
# decomposition of the RCG-04E wide conical-spiral state -- ``bs1``/``bs2``
# are a fixed-physical-partition spatial-refinement pair (block_size_x 1 vs
# 2, same 48-fine/32-interface/112-coarse-atom split), ``bs1_shifted`` moves
# the FINE seed half a period away at the same block resolution. See
# run_moving_static_mixed.py, static_topology_oracle.py, and each fixture's
# README.md.
MOVING_STATIC_MIXED_BS1_CASE = "moving_static_mixed_bs1"
MOVING_STATIC_MIXED_BS2_CASE = "moving_static_mixed_bs2"
MOVING_STATIC_MIXED_BS1_SHIFTED_CASE = "moving_static_mixed_bs1_shifted"
MOVING_STATIC_MIXED_CASES = (
    MOVING_STATIC_MIXED_BS1_CASE, MOVING_STATIC_MIXED_BS2_CASE,
    MOVING_STATIC_MIXED_BS1_SHIFTED_CASE,
)

# RCG-04G E2E-MOVING-ADAPTIVE: a deterministic periodic domain-wall pair
# (RCG-04B's domain_wall_pair_state, Initmag=4) under ADAPTIVE mask mode,
# demonstrating accepted resolution transitions correlated with genuine
# wall motion and block-boundary crossing, plus the RCG-03 polarization
# safety contract. ``moving_wall_feature_off`` is the plain (AdaptiveCG
# disabled) physics reference sharing a byte-identical restart_seed.out.
# See run_moving_adaptive_wall.py and each fixture's README.md.
MOVING_WALL_FEATURE_OFF_CASE = "moving_wall_feature_off"
MOVING_WALL_ADAPTIVE_CASE = "moving_wall_adaptive"

# RCG-04H DMI-HYBRID-CROSSING: RCG-04B's deterministic +q/-q chiral partner
# pair (Initmag=8) on the RCG-04F fine/interface/coarse partition, with DMI
# (``dmfile_chiral``) and uniaxial anisotropy (``kfile_cg_x``) both enabled,
# asserting signed chirality, +q/-q DMI energy ordering, named DMI/
# anisotropy energy series, a spatial-refinement pair, and DMI-specific
# interface ownership -- plus a sign-reversed-DMI negative control
# (``dmfile_chiral_reversed``, a tracked exact ``D -> -D`` negation, not a
# temporary source mutation). See run_moving_dmi_chiral.py and each
# fixture's README.md.
MOVING_DMI_CHIRAL_ALL_FINE_PLUS_CASE = "moving_dmi_chiral_all_fine_plus"
MOVING_DMI_CHIRAL_BS1_PLUS_CASE = "moving_dmi_chiral_bs1_plus"
MOVING_DMI_CHIRAL_BS1_MINUS_CASE = "moving_dmi_chiral_bs1_minus"
MOVING_DMI_CHIRAL_BS2_PLUS_CASE = "moving_dmi_chiral_bs2_plus"
MOVING_DMI_CHIRAL_BS1_PLUS_REVERSED_CASE = "moving_dmi_chiral_bs1_plus_reversed"
MOVING_DMI_CHIRAL_BS1_MINUS_REVERSED_CASE = "moving_dmi_chiral_bs1_minus_reversed"
MOVING_DMI_CHIRAL_CASES = (
    MOVING_DMI_CHIRAL_ALL_FINE_PLUS_CASE, MOVING_DMI_CHIRAL_BS1_PLUS_CASE,
    MOVING_DMI_CHIRAL_BS1_MINUS_CASE, MOVING_DMI_CHIRAL_BS2_PLUS_CASE,
    MOVING_DMI_CHIRAL_BS1_PLUS_REVERSED_CASE, MOVING_DMI_CHIRAL_BS1_MINUS_REVERSED_CASE,
)

# RCG-05C OWNERSHIP-ANISO-BUFFER: a genuinely anisotropic-block-shape
# (`buffer_width_blocks=(2,1,1)`), `cg_mask_mode ADAPTIVE` fixture that
# engages the real GPU adaptive-transition runtime, used by
# `run_ownership_map_comparator.py` to demonstrate the GPU buffer-width
# scalarization defect through the real executable. See
# `e2e/ownership_aniso_buffer/README.md`.
OWNERSHIP_ANISO_BUFFER_CASE = "ownership_aniso_buffer"


def all_e2e_cases() -> tuple[str, ...]:
    """Return every e2e directory referenced by either executable harness."""
    names = {
        FEATURE_OFF_CASE, STATIC_ALL_FINE_CASE, STATIC_ALL_COARSE_CASE,
        STATIC_MIXED_CASE, ADAPTIVE_MIXED_CASE, POLARIZATION_GATE_CASE,
        POLARIZATION_GATE_GPU_CASE, *CPU_REJECTION_CASES,
        *CPU_ANISOTROPY_CASES, CLUSTER_ANISOTROPY_REJECTION_CASE, *CPU_PARITY_CASES.values(),
        CPU_INITIAL_PHASE_SD_CASE, *[case for case, _ in CPU_INITIAL_PHASE_CASES.values()],
        SPIN_SPIRAL_CASE, RESTART_INITMAG_CASE, *GPU_EXECUTABLE_CASES, *GPU_INITIAL_PHASE_CASES,
        *GPU_PARITY_CASES.values(), GPU_DMI_CASE, GPU_FFT_CASE,
        MOVING_FEATURE_OFF_CASE, MOVING_ALL_FINE_CASE,
        MOVING_FEATURE_OFF_WIDE_CASE, MOVING_ALL_FINE_WIDE_CASE, *MOVING_ALL_COARSE_CASES,
        *MOVING_STATIC_MIXED_CASES,
        MOVING_WALL_FEATURE_OFF_CASE, MOVING_WALL_ADAPTIVE_CASE,
        *MOVING_DMI_CHIRAL_CASES,
        OWNERSHIP_ANISO_BUFFER_CASE,
    }
    return tuple(sorted(names))
