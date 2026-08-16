#pragma once

#include "tensor.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "gpuStructures.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpu_wrappers.h"
using ParallelizationHelper = GpuParallelizationHelper;

// CGP-03B: an adaptive-CG-owned sibling of GpuMomentUpdater.
//
// GpuMomentUpdater itself (source/gpu_files/gpuMomentUpdater.{hpp,cpp}) is
// shared, non-adaptive infrastructure that every SD path -- adaptive or not
// -- calls unconditionally, once per step, over every atom. CGP-03's
// evidence doc (docs/CGP-03_LAZY_MATERIALIZATION_EVIDENCE.md, section 2)
// found that this is the hard blocker standing between CGP-03's "zero
// full-lattice reconstruction on an ordinary step" target and reality: even
// once atom-direction commit becomes lazy, GpuMomentUpdater's mmom2/emomM
// kernels still sweep every atom, unconditionally, regardless of adaptive
// state.
//
// This class gives adaptive dynamics its own moment-update entry point
// instead of editing GpuMomentUpdater. GpuMomentUpdater is not modified by
// this task in any way -- see gpuMomentUpdaterSourceHash in the test suite
// for the negative control that proves it.
//
// See docs/CGP-03B_ADAPTIVE_MOMENT_UPDATER_EVIDENCE.md for the full
// audit, contract, and evidence.
class GpuAdaptiveMomentUpdater {
private:
   // Moments to update -- identical binding to GpuMomentUpdater's, since
   // both update the one shared deviceLattice storage (see the contract's
   // "same buffers, different launch geometry" note in the evidence doc).
   GpuTensor<real, 2>& mmom;
   GpuTensor<real, 2>& mmom0;
   GpuTensor<real, 2>& mmom2;
   GpuTensor<real, 3>& emom;
   GpuTensor<real, 3>& emom2;
   GpuTensor<real, 3>& emomM;
   GpuTensor<real, 2>& mmomi;

   // Parameters
   int mompar;
   char initexc;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   ParallelizationHelper& parallel;

public:
   // Parallelization classes.  Deliberately new types, not aliases of
   // GpuMomentUpdater::Mompar1/Mompar2/Copy1/Copy2 -- those are private
   // members of GpuMomentUpdater, and reaching into them (via friend or
   // otherwise) would mean touching gpuMomentUpdater.hpp, which this task's
   // design requirement rules out. The four per-atom formulas below are the
   // same math, intentionally duplicated in one place rather than shared,
   // and are guarded against drift by testAdaptiveMomentUpdaterFullParity
   // (bitwise cross-check against the real GpuMomentUpdater on every call).
   class Mompar1;
   class Mompar2;
   class Copy1;
   class Copy2;

   // Constructor -- same shape as GpuMomentUpdater's.
   GpuAdaptiveMomentUpdater(deviceLattice& gpuLattice, int mompar, char initexc);

   // Full-lattice update. Every kernel launch geometry, formula, and buffer
   // touched is identical to GpuMomentUpdater::update(). CGP-03B could not
   // safely wire updateActiveOnly() into production because
   // GpuCorrelations::measure()'s firing predicate (source/gpu_files/
   // correlations/gpuCorrelations.cpp) mutates its own bookkeeping
   // (t_cur, curstep-derived state) as a side effect, so it cannot be
   // "peeked" one step ahead without editing shared, non-adaptive code.
   // CGP-03C's production wiring (gpuSDSimulation.cpp, via
   // adaptiveMomentUpdateNeedsFullLattice() below) still calls this
   // wherever that look-ahead determines a full update is needed -- see
   // docs/CGP-03C_MOMENT_UPDATE_LOOKAHEAD_EVIDENCE.md.
   void updateFull();

   // Touched-atom-only update: the Mompar1/Mompar2 (mmom2 rebuild) and
   // Copy1/Copy2 (mmomi/emomM rebuild) kernels run only over
   // activeAtoms/activeCount -- the same compact, one-based, per-ensemble
   // list already driving GpuDepondtIntegrator::evolveFirst/evolveSecond for
   // adaptive steps (GpuAdaptiveRuntime::activeAtomList()/activeAtomCount())
   // -- rather than every atom in the system.
   //
   // emom<->emom2 and mmom<->mmom2 stay whole-tensor pointer swaps (O(1),
   // GpuTensor::swap, unchanged cost and unchanged for every atom): the
   // atom-direction contract (CGP-03) is untouched by this task, and a
   // partial-content swap would be both unnecessary (the swap itself is
   // free) and unsafe (it would leave emom/emom2 inconsistent between
   // touched and untouched atoms across the very next materialization).
   //
   // Atoms not in the active list keep whatever mmom/mmomi/emomM value they
   // held before this call -- by construction, not by omission: neither
   // buffer slot for that atom is ever written, so the same numeric value
   // simply keeps changing which tensor (mmom vs mmom2) currently names it.
   //
   // CGP-03C: this is now the production default for adaptive runs on
   // steps where adaptiveMomentUpdateNeedsFullLattice() below returns
   // false; the caller must pass a *freshly re-fetched*
   // GpuAdaptiveRuntime::activeAtomList()/activeAtomCount() (not a snapshot
   // taken before a transition may have been published this same step).
   void updateActiveOnly(const int* activeAtoms, std::size_t activeCount);
};

// Forward-declared rather than #include "gpuMeasurement.hpp" here, to keep
// this header's dependency footprint the same as before CGP-03C for every
// consumer that does not need the look-ahead.
class GpuMeasurement;

// CGP-03C: pure, side-effect-free query for whether the *next* loop
// iteration's measurement/output traffic will need gpuLattice.mmom/mmomi/
// emomM fully fresh for every atom -- i.e. whether the current step's
// adaptive moment update must be GpuAdaptiveMomentUpdater::updateFull()
// rather than updateActiveOnly(). Deliberately does not consider
// GpuCorrelations/do_sc at all -- see
// docs/CGP-03C_MOMENT_UPDATE_LOOKAHEAD_EVIDENCE.md for why that is an
// accepted, deliberate scope decision (Anders: correlations are not
// guaranteed correct under adaptive CG), not an oversight.
//
// gpuMeasurement may be null (MeasurementFactory chose FortranMeasurement,
// i.e. do_cuda_measurements != 'Y'); the Fortran do_measurements oracle
// alone then matches that backend's own existing freshness guarantee, so no
// special-casing is needed for it.
//
// nstep/rstep/currentStep/lastStep are plain values (not a GpuSimulation&)
// so this function has no dependency on GpuSimulation's internals and can
// be unit-tested standalone.
bool adaptiveMomentUpdateNeedsFullLattice(GpuMeasurement* gpuMeasurement,
                                          std::size_t currentStep, std::size_t lastStep,
                                          std::size_t nstep, std::size_t rstep);
