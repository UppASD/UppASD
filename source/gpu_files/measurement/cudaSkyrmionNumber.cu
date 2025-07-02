#include "cudaSkyrmionNumber.cuh"
#include "fortranData.hpp"
#include "cudaParallelizationHelper.hpp"


__global__ void grad_moments_kernel(const CudaTensor<real, 3> emomM,
                                    const CudaTensor<real, 3> dxyz_vec,
                                    const CudaTensor<int, 2> dxyz_atom,
                                    const CudaTensor<int, 1> dxyz_list,
                                    CudaTensor<real, 4> grad_mom);

__global__ void pontryagin_no_kernel(const CudaTensor<real, 3> emomM,
                                     const CudaTensor<real, 4> grad_mom,
                                     real* pontryagin_no_out);

__global__ void avrg_skyno_kernel();


SkyrmionNumber::SkyrmionNumber(const CudaTensor<real, 3>& emomM, const CudaTensor<real, 3>& emom)
: emomM(emomM)
, emom(emom)
, skyno_step(*FortranData::skyno_step)
, buffer_size(*FortranData::skyno_buff)
{
    const uint N = emomM.extent(1);
    const uint M = emomM.extent(2);

    Tensor<real, 3> dxyz_vec_fortran(FortranData::dxyz_vec, 3, 26, N);
    dxyz_vec.Allocate(3, 26, N); // why 26?
    dxyz_vec.copy_sync(dxyz_vec_fortran);

    Tensor<int, 2> dxyz_atom_fortran(FortranData::dxyz_atom, 26, N);
    dxyz_atom.Allocate(26, N);
    dxyz_atom.copy_sync(dxyz_atom_fortran);

    Tensor<int, 1> dxyz_list_fortran(FortranData::dxyz_list, N);
    dxyz_list.Allocate(N);
    dxyz_list.copy_sync(dxyz_list_fortran);

//    std::cout << "dxyz_vec\n";
//    for (uint i = 0; i < N; ++i)
//    {
//        std::cout << i << ": ";
//        for (uint k = 0; k < 26; ++k)
//        {
//            std::cout << "(" << dxyz_vec_fortran(0, k, i) << " "
//                      << dxyz_vec_fortran(1, k, i) << " "
//                      << dxyz_vec_fortran(2, k, i) << "), ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "\ndxyz_atom\n";
//    for (uint i = 0; i < N; ++i)
//    {
//        std::cout << i << ": ";
//        for (uint k = 0; k < 26; ++k)
//        {
//            std::cout << dxyz_atom_fortran(k, i) << " ";
//        }
//        std::cout << "\n";
//    }
//
//    std::cout << "\ndxyz_list\n";
//    for (uint i = 0; i < N; ++i)
//    {
//        std::cout << i << ": " << dxyz_list_fortran(i) << "\n";
//    }

    grad_mom.Allocate(3, 3, N, M);
    grad_mom.zeros();

    skynob.Allocate(buffer_size);
    skynob.zeros();

    sk_avrg.Allocate(buffer_size);
    sk_avrg.zeros();

    sk_var.Allocate(buffer_size);
    sk_var.zeros();

    // TODO: this should not be allocated, but hooked up to fortran buffer pointer
    indxb_skyno.AllocateHost(buffer_size);
}


SkyrmionNumber::~SkyrmionNumber()
{
    dxyz_vec.Free();
    dxyz_atom.Free();
    dxyz_list.Free();
    grad_mom.Free();
    skynob.Free();
    sk_avrg.Free();
    sk_var.Free();
    indxb_skyno.FreeHost();
}


void SkyrmionNumber::measure(std::size_t mstep)
{
    --mstep;

    if (mstep % skyno_step != 0)
        return;

    // std::cout << "[SkyrmionNumber::measure] mstep = " << mstep << ", ";

    const uint N = emomM.extent(1);
    const uint M = emomM.extent(2);

    // this seems to be done on every time step in fortran
    cudaStream_t workStream = CudaParallelizationHelper::def.getWorkStream();
    dim3 threads = {512};
    dim3 blocks = {
            (N + threads.x - 1) / threads.x,
            (M + threads.y - 1) / threads.y
    };

    // TODO are we sure this should be emom, and not emomM?
    grad_moments_kernel<<<blocks, threads, 0, workStream>>>(emom, dxyz_vec, dxyz_atom, dxyz_list, grad_mom);
    pontryagin_no_kernel<<<blocks, threads, 0, workStream>>>(emomM, grad_mom, skynob.data() + buffer_count);
    indxb_skyno(buffer_count++) = static_cast<uint>(mstep);

    // cudaDeviceSynchronize();

    if (buffer_count >= buffer_size)
    {
        // TODO: copy to fortran
        buffer_count = 0;
        skynob.zeros();
    }
}


void SkyrmionNumber::flushMeasurements(std::size_t mstep)
{

}


__global__ void grad_moments_kernel(const CudaTensor<real, 3> emomM,
                                    const CudaTensor<real, 3> dxyz_vec,
                                    const CudaTensor<int, 2> dxyz_atom,
                                    const CudaTensor<int, 1> dxyz_list,
                                    CudaTensor<real, 4> grad_mom)
{
    const uint N = emomM.extent(1);
    const uint M = emomM.extent(2);
    const uint iatom = blockDim.x * blockIdx.x + threadIdx.x;
    const uint kk = blockDim.y * blockIdx.y + threadIdx.y;

    if (iatom >= N || kk >= M)
        return;

    assert(dxyz_list(iatom) < N);

    for (uint jneigh = 0; jneigh < dxyz_list(iatom); ++jneigh)
    {
        assert(jneigh < 26);
        const uint jatom = dxyz_atom(jneigh, iatom) - 1; // needs -1 here since it gives the index of a neighboring atom

        assert(jatom < N);
        const real d_mom[3] = {
            emomM(0, jatom, kk) - emomM(0, iatom, kk),
            emomM(1, jatom, kk) - emomM(1, iatom, kk),
            emomM(2, jatom, kk) - emomM(2, iatom, kk)
        };

        const real dv[3] = { // dv = {dx, dy, dz}
            dxyz_vec(0, jneigh, iatom),
            dxyz_vec(1, jneigh, iatom),
            dxyz_vec(2, jneigh, iatom)
        };

        for (uint coord = 0; coord < 3; ++coord)
        {
            if (abs( dv[coord] ) > 1e-7)
            {
                grad_mom(0, coord, iatom, kk) += d_mom[0] / dv[coord];
                grad_mom(1, coord, iatom, kk) += d_mom[1] / dv[coord];
                grad_mom(2, coord, iatom, kk) += d_mom[2] / dv[coord];
            }
        }
    }

    for (uint coord = 0; coord < 3; ++coord)
    {
        grad_mom(coord, 0, iatom, kk) /= dxyz_list(iatom);
        grad_mom(coord, 1, iatom, kk) /= dxyz_list(iatom);
        grad_mom(coord, 2, iatom, kk) /= dxyz_list(iatom);
    }
}


__global__ void pontryagin_no_kernel(const CudaTensor<real, 3> emomM,
                                     const CudaTensor<real, 4> grad_mom,
                                     real* pontryagin_no_out)
{
    const uint N = emomM.extent(1);
    const uint M = emomM.extent(2);

    const uint iatom = blockDim.x * blockIdx.x + threadIdx.x;
    const uint k = blockDim.y * blockIdx.y + threadIdx.y;

    if (iatom >= N || k >= M)
        return;

    const real cvec_x = grad_mom(1,0,iatom,k) * grad_mom(2,1,iatom,k)
                      - grad_mom(2,0,iatom,k) * grad_mom(1,1,iatom,k);

    const real cvec_y = grad_mom(2,0,iatom,k) * grad_mom(0,1,iatom,k)
                      - grad_mom(0,0,iatom,k) * grad_mom(2,1,iatom,k);

    const real cvec_z = grad_mom(0,0,iatom,k) * grad_mom(1,1,iatom,k)
                      - grad_mom(1,0,iatom,k) * grad_mom(0,1,iatom,k);


    const real partial_sum = emomM(0,iatom,k) * cvec_x
                           + emomM(1,iatom,k) * cvec_y
                           + emomM(2,iatom,k) * cvec_z;

    atomicAdd(pontryagin_no_out, partial_sum);

    __syncthreads();

    if (iatom == 0 && k == 0)
    {
        *pontryagin_no_out /= (M_PI * M);
    }
}


__global__ void avrg_skyno_kernel(CudaTensor<real, 1> skynob,
                                  CudaTensor<real, 1> sk_avrg,
                                  CudaTensor<real, 1> sk_var)
{
    const uint buffer_count = skynob.extent(0);
    for (uint k = 0; k < buffer_count; ++k)
    {

    }
}














