#include "depondt_kernels.hpp"


__global__ void InitializeDeviceRngStates(GpuRngState* states, unsigned long long seed, unsigned int count) {

   const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
   if (id < count) {
      GPU_RAND_INIT(seed, static_cast<unsigned long long>(id), 0ULL, &states[id]);
   }
}

__global__ void InitializeThermalSigmaFactor(GpuTensor<real, 1> sigmaFactor, unsigned int count, real dp_factor) {

   const unsigned int site = blockIdx.x * blockDim.x + threadIdx.x;
   if (site < count) {
      sigmaFactor[site] = sqrt(dp_factor * sigmaFactor[site]);
   }
}

__device__ inline void rotateDepondtVector(real mx, real my, real mz, real bx, real by, real bz, real timestep, real gamma, 
                                        real damping, real& outx, real& outy, real& outz){
    const real norm = sqrt(bx * bx + by * by + bz * bz);

    if (norm == real(0)) {
        outx = mx;
        outy = my;
        outz = mz;
        return;
    }
    bx /= norm;
    by /= norm;
    bz /= norm;

    real angle = norm * timestep * gamma;
    angle /= real(1.0) + damping * damping;

    real cosv, sinv;
    sincos(angle, &sinv, &cosv);
    const real u = real(1.0) - cosv;

    outx = mx * (bx * bx * u + cosv)
         + my * (bx * by * u - bz * sinv)
         + mz * (bx * bz * u + by * sinv);

    outy = mx * (bx * by * u + bz * sinv)
         + my * (by * by * u + cosv)
         + mz * (by * bz * u - bx * sinv);

    outz = mx * (bx * bz * u - by * sinv)
         + my * (bz * by * u + bx * sinv)
         + mz * (bz * bz * u + cosv);
}

__device__ inline void BuildEffectiveField(real& ex, real& ey, real& ez, const real mx, const real my, const real mz,
                                           const real lx, const real ly, const real lz, const real damping){

    ex = lx + damping * (my * lz - mz * ly);
    ey = ly + damping * (mz * lx - mx * lz);
    ez = lz + damping * (mx * ly - my * lx);
}

__device__ inline void RandomField(real& x, real& y, real& z, GpuRngState* d_state, unsigned int atom, unsigned int site, unsigned int ensemble,
                                 GpuTensor<real, 1> sigmaFactor, GpuTensor<real, 2> mmom){
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

    const real sigma = sigmaFactor(site) * rsqrt(mmom(site, ensemble));

    x *= sigma;
    y *= sigma;
    z *= sigma;
}

__global__ void DepondtEvolveFirst(GpuTensor<real, 3> thermfield, const GpuTensor<real, 1> sigmaFactor, GpuTensor<GPU_RAND_STATE, 1> d_state,
                                  const GpuTensor<real, 3> beff, GpuTensor<real, 3> b2eff, GpuTensor<real, 3> emom, GpuTensor<real, 3> emom2, GpuTensor<real, 2> mmom, 
                                  GpuTensor<real, 3> emomM, const GpuTensor<real, 3> btorque, int N, int M,  real timestep, real gamma, real damping, bool applyTorque){


    const unsigned int site = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int ensemble = blockIdx.y;

    if (site >= N || ensemble >= M)
        return;

    const unsigned int atom = ensemble * N + site;

    real tx, ty, tz;
    real ex, ey, ez;
    real rx, ry, rz;

    const real mx = emom(0, site, ensemble);
    const real my = emom(1, site, ensemble);
    const real mz = emom(2, site, ensemble);

    const real bx = beff(0, site, ensemble);
    const real by = beff(1, site, ensemble);
    const real bz = beff(2, site, ensemble);

    RandomField(tx, ty, tz, d_state.data(), atom, site, ensemble, sigmaFactor, mmom);

    thermfield(0, site, ensemble) = tx;
    thermfield(1, site, ensemble) = ty;
    thermfield(2, site, ensemble) = tz;

    const real lx = bx + tx;
    const real ly = by + ty;
    const real lz = bz + tz;

    BuildEffectiveField(ex, ey, ez, mx, my, mz, lx, ly, lz, damping);

    if (applyTorque) {
        ex += btorque(0, site, ensemble);
        ey += btorque(1, site, ensemble);
        ez += btorque(2, site, ensemble);
    }

    b2eff(0, site, ensemble) = ex;
    b2eff(1, site, ensemble) = ey;
    b2eff(2, site, ensemble) = ez;

    rotateDepondtVector(mx, my, mz, ex, ey, ez, timestep, gamma, damping, rx, ry, rz);

    emom2(0, site, ensemble) = mx;
    emom2(1, site, ensemble) = my;
    emom2(2, site, ensemble) = mz;

    emom(0, site, ensemble) = rx;
    emom(1, site, ensemble) = ry;
    emom(2, site, ensemble) = rz;

    const real moment = mmom(site, ensemble);

    emomM(0, site, ensemble) = rx * moment;
    emomM(1, site, ensemble) = ry * moment;
    emomM(2, site, ensemble) = rz * moment;

   
};

__global__ void DepondtEvolveSecond(const GpuTensor<real, 3> thermfield,const GpuTensor<real, 3> beff, const GpuTensor<real, 3> b2eff,
                                    GpuTensor<real, 3> emom, GpuTensor<real, 3> emom2, const GpuTensor<real, 3> btorque,
                                    int N, int M, real timestep, real gamma, real damping, bool applyTorque){
    const unsigned int site = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int ensemble = blockIdx.y;

    if (site >= N || ensemble >= M)
        return;

    real ex, ey, ez;
    real rx, ry, rz;

    const real mx = emom(0, site, ensemble);
    const real my = emom(1, site, ensemble);
    const real mz = emom(2, site, ensemble);

    const real ox = emom2(0, site, ensemble);
    const real oy = emom2(1, site, ensemble);
    const real oz = emom2(2, site, ensemble);

    const real bx = beff(0, site, ensemble);
    const real by = beff(1, site, ensemble);
    const real bz = beff(2, site, ensemble);

    const real lx = bx + thermfield(0, site, ensemble);
    const real ly = by + thermfield(1, site, ensemble);
    const real lz = bz + thermfield(2, site, ensemble);

    BuildEffectiveField(ex, ey, ez,mx, my, mz, lx, ly, lz,damping);

    if (applyTorque) {
        ex += btorque(0, site, ensemble);
        ey += btorque(1, site, ensemble);
        ez += btorque(2, site, ensemble);
    }

    ex = real(0.5) * (ex + b2eff(0, site, ensemble));
    ey = real(0.5) * (ey + b2eff(1, site, ensemble));
    ez = real(0.5) * (ez + b2eff(2, site, ensemble));

    rotateDepondtVector(ox, oy, oz,ex, ey, ez,timestep, gamma, damping, rx, ry, rz);

    emom(0, site, ensemble) = ox;
    emom(1, site, ensemble) = oy;
    emom(2, site, ensemble) = oz;

    emom2(0, site, ensemble) = rx;
    emom2(1, site, ensemble) = ry;
    emom2(2, site, ensemble) = rz;
}

