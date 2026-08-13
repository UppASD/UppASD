#pragma once

// RCG-06D (F-20): shared tuple-seeded MINSTD reconstruction RNG.
//
// Constrained-cone reconstruction seeds a per-(block,channel,ensemble,epoch)
// stream deterministically from adaptive_reconstruction_configuration_type's
// global_seed, so a replay is reproducible without storing per-atom random
// state. This slice found and fixed two independent CPU/GPU divergences in
// this device-side generator, both invisible to every tracked fixture
// because they all use RECONSTRUCTION_ALIGNED (or CONSTRAINED_CONE at
// cone_angle_rad=0), where the draw sequence is unconditionally nullified
// (multiplied by an amplitude that bisection forces to exactly zero) before
// it can affect atom_direction -- see the RCG-06D evidence doc:
//
// 1. Multiplier mismatch: this generator used 16807 (the original
//    Lehmer/"MINSTD" constant) while the CPU implementation
//    (adaptivehybridsolver.f90's next_uniform) used 48271 (the revised
//    Park-Miller "minimal standard" constant), both against the same
//    modulus 2^31-1. Fixed to 48271 unconditionally, matching CPU.
// 2. Seed-formula epoch offset: block/channel/ensemble are 0-based on the
//    GPU call site (gpuAdaptiveRuntime.cpp's loop indices) vs 1-based in
//    Fortran, so the "+1" applied to each of those three here is a genuine,
//    necessary index-convention conversion -- confirmed correct. epoch has
//    no such convention difference: both sides pass the same
//    zero-based, monotonically-incrementing transition-epoch counter
//    (blockselector.f90's transition_epoch / gpuAdaptiveRuntime.hpp's
//    transitionEpoch, both initialized to 0 and incremented by exactly 1 on
//    each accepted transition). Unconditionally adding 1 to epoch here too
//    was therefore not a convention conversion but a plain off-by-one:
//    GPU's epoch=0 silently matched what CPU's formula calls epoch=1,
//    permanently offsetting every GPU reconstruction stream by one epoch
//    relative to the CPU stream it is meant to replay. Fixed by dropping
//    the "+1" on epoch only.
//
// Shared between gpuAdaptiveRuntime.cpp (the production kernels) and
// test_reconstruction_rng_gpu_parity.cpp (the standalone RCG-06D CPU/GPU
// parity fixture), so the fixture measures the literal generator the
// production kernels use, not a reimplementation that could silently drift.

#include "real_type.h"

constexpr std::uint64_t adaptiveReconstructionModulus = 2147483647ULL;

__device__ inline std::uint64_t adaptiveTupleSeed(std::uint64_t globalSeed, std::size_t block,
                                                  std::size_t channel, std::size_t ensemble,
                                                  unsigned int epoch) {
   std::uint64_t seed = globalSeed % adaptiveReconstructionModulus;
   seed = (seed + 104729ULL * ((block + 1) % adaptiveReconstructionModulus)) %
      adaptiveReconstructionModulus;
   seed = (seed + 130363ULL * ((channel + 1) % adaptiveReconstructionModulus)) %
      adaptiveReconstructionModulus;
   seed = (seed + 155921ULL * ((ensemble + 1) % adaptiveReconstructionModulus)) %
      adaptiveReconstructionModulus;
   seed = (seed + 196613ULL * (static_cast<std::uint64_t>(epoch) %
                              adaptiveReconstructionModulus)) % adaptiveReconstructionModulus;
   if(seed == 0) seed = adaptiveReconstructionModulus - 1;
   return seed;
}

__device__ inline real adaptiveNextUniform(std::uint64_t& state) {
   state = (48271ULL * state) % adaptiveReconstructionModulus;
   return static_cast<real>(state) / static_cast<real>(adaptiveReconstructionModulus);
}
