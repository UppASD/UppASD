#include "depondt_kernals.hpp"


__global__ void InitializeDeviceRngStates(GpuRngState* states, unsigned long long seed, unsigned int count) {

   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if (id < count) {
      GPU_RAND_INIT(seed, static_cast<unsigned long long>(id), 0ULL, &states[id]);
   }
}

__global__ void InitializeThermalSigmaFactor(GpuTensor<real, 1> sigmaFactor, const GpuTensor<real, 1> temperature,
                                             unsigned int count, real damping) {

   const unsigned int site = blockIdx.x * blockDim.x + threadIdx.x;
   if (site < count) {
      sigmaFactor[site] = sqrt(damping * temperature[site]);
   }
}

__device__ inline void rotateDepondtVector(){}

__device__ inline void BuildEffectiveField(real& bx, real& by, real& bz, const real mx, const real my, const real mz,
                                           const real lx, const real ly, const real lz, const real damping){

    bx = lx + damping * (my * lz - mz * ly);
    by = ly + damping * (mz * lx - mx * lz);
    bz = lz + damping * (mx * ly - my * lx);
}

__device__ inline void RandomField(real& x, real& y, real& z, GPU_RAND_STATE* d_state, unsigned int atom, unsigned int site, unsigned int ensemble,
                                 GpuTensor<real, 1> sigmaFactor, GpuTensor<real, 2> mmom)
{
    GPU_RAND_STATE state = d_state[atom];

#if defined(SINGLE_PREC)
    x = GPU_NORMAL(&state);
    y = GPU_NORMAL(&state);
    z = GPU_NORMAL(&state);
#else
    x = GPU_NORMAL_DOUBLE(&state);
    y = GPU_NORMAL_DOUBLE(&state);
    z = GPU_NORMAL_DOUBLE(&state);
#endif

    d_state[atom] = state;
.
    const real sigma = sigmaFactor(site) * rsqrt(mmom(site, ensemble));

    x *= sigma;
    y *= sigma;
    z *= sigma;
}

__global__ void DepondtEvolveFirst(GpuTensor<real, 3> thermfield, const GpuTensor<real, 1> sigmaFactor, GpuTensor<GPU_RAND_STATE, 1> d_state,
                                  const GpuTensor<real, 3> beff, GpuTensor<real, 3> emom, GpuTensor<real, 3> emom2, GpuTensor<real, 2> mmom, 
                                  GpuTensor<real, 3> emomM, const GpuTensor<real, 3> btorque, int N, int M, bool applyTorque){

const unsigned int site = blockIdx.x * blockDim.x + threadIdx.x;
const unsigned int ensemble = blockIdx.y;

if (site >= N) return;

const unsigned int atom = ensemble * N + site;

real tx, ty, tz, bx, by, bz, lx, ly, lz;
real mx, my, mz;

const mx = emom(0, site, ensemble);
const my = emom(1, site, ensemble);
const mz = emom(2, site, ensemble);

const real bx = beff(0, site, ensemble);
const real by = beff(1, site, ensemble);
const real bz = beff(2, site, ensemble);

if (applyTorque) {
    bx += btorque(0, site, ensemble);
    by += btorque(1, site, ensemble);
    bz += btorque(2, site, ensemble);
}

BuildRandomField(tx, ty, tz, d_state.data(), atom, site, ensemble, sigmaFactor.data(), mmom.data());

thermfield(0, site, ensemble) = tx;
thermfield(1, site, ensemble) = ty;
thermfield(2, site, ensemble) = tz;

lx = bx + tx;
ly = by + ty;
lz = bz + tz;

BuildEffectiveField(ex, ey, ez, mx, my, mz, lx, ly, lz, damping);
   


   
};

__global__ void DepondtEvolveSecond(){}; 

