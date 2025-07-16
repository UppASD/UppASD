#include "cudaMeasurement.cuh"

#include <cuda_runtime.h>

#include "c_helper.h"
#include "stopwatchPool.hpp"
#include "fortranData.hpp"
#include "cudaParallelizationHelper.hpp"


CudaMeasurement::CudaMeasurement(const CudaTensor<real, 3>& emomM,
                                 const CudaTensor<real, 3>& emom,
                                 const CudaTensor<real, 2>& mmom,
                                 bool alwaysCopy)
: emomM(emomM)
, emom(emom)
, mmom(mmom)
, stopwatch(GlobalStopwatchPool::get("Cuda measurement"))
{
    if (*FortranData::do_avrg == 'Y')
    {
        assert(*FortranData::avrg_step > 0 && *FortranData::avrg_buff > 0);

        mavg_buff_gpu.Allocate(*FortranData::avrg_buff);
        mavg_buff_gpu.zeros();

        mavg_buff_cpu.AllocateHost(*FortranData::avrg_buff);
        mavg_buff_cpu.zeros();

        mavg_iter.AllocateHost(*FortranData::avrg_buff);
        mavg_iter.zeros();

        std::cout << "AverageMagnetization observable added" << std::endl;
    }

    if (*FortranData::do_cumu == 'Y')
    {
        assert(*FortranData::cumu_step > 0 && *FortranData::cumu_buff > 0);

        cumu_buff_gpu.Allocate(1);
        cumu_buff_gpu.zeros();

        cumu_buff_cpu.AllocateHost(1);
        cumu_buff_cpu.zeros();

        std::cout << "BinderCumulant observable added" << std::endl;
    }

    if (*FortranData::do_skyno == 'Y')
    {
        assert(*FortranData::skyno_step > 0 && *FortranData::skyno_buff > 0);

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

        grad_mom.Allocate(3, 3, N, M);
        grad_mom.zeros();

        skyno_buff_gpu.Allocate(3, *FortranData::skyno_buff);
        skyno_buff_gpu.zeros();

        skyno_buff_cpu.AllocateHost(3, *FortranData::skyno_buff);
        skyno_buff_cpu.zeros();

        skyno_iter.AllocateHost(*FortranData::skyno_buff);
        skyno_iter.zeros();

        std::cout << "SkyrmionNumber observable added" << std::endl;
    }
}


CudaMeasurement::~CudaMeasurement()
{
    if (*FortranData::do_avrg == 'Y')
    {
        mavg_buff_gpu.Free();
        mavg_buff_cpu.FreeHost();
        mavg_iter.FreeHost();
    }


    if (*FortranData::do_cumu == 'Y')
    {
        cumu_buff_gpu.Free();
        cumu_buff_cpu.FreeHost();
    }

    if (*FortranData::do_skyno == 'Y')
    {
        dxyz_vec.Free();
        dxyz_atom.Free();
        dxyz_list.Free();
        grad_mom.Free();
        skyno_buff_gpu.Free();
        skyno_buff_cpu.FreeHost();
        skyno_iter.FreeHost();
    }
}


void CudaMeasurement::measure(std::size_t mstep)
{
    --mstep; // this is because the simulation loop begins at 1 because of Fortran indexing

    if (*FortranData::do_avrg == 'Y')
        measureAverageMagnetization(mstep);

    if (*FortranData::do_cumu == 'Y')
        measureBinderCumulant(mstep);

    if (*FortranData::do_skyno == 'Y')
        measureSkyrmionNumber(mstep);
}


void CudaMeasurement::flushMeasurements(std::size_t mstep)
{
 // TODO: perhaps if alwaysCopy is false we could queue up measurements from multiple
 // mstep's and then only sync if buffer is full or if flushMeasurements is called
}



__global__ void naiveAverageMagnetization_kernel(const CudaTensor<real, 3> emomM,
                                                 CudaMeasurement::AverageMagnetizationData* d)
{
    const uint k = blockIdx.x * blockDim.x + threadIdx.x;
    const uint atoms = emomM.extent(1);
    const uint ensembles = emomM.extent(2);

    if (k >= ensembles)
        return;

    real m[3] = {0};

    for (uint i = 0; i < atoms; ++i)
    {
        m[0] += emomM(0, i, k);
        m[1] += emomM(1, i, k);
        m[2] += emomM(2, i, k);
    }

    m[0] /= atoms;
    m[1] /= atoms;
    m[2] /= atoms;

    const real m_norm = norm(3, m);
    const real avrgms = pow(m_norm, 2);

    atomicAdd(&(d->m_x), m[0]);
    atomicAdd(&(d->m_y), m[1]);
    atomicAdd(&(d->m_z), m[2]);
    atomicAdd(&(d->m), m_norm);
    atomicAdd(&(d->m_stdv), avrgms);

    __syncthreads();

    if (k == 0)
    {
        d->m_x /= ensembles;
        d->m_y /= ensembles;
        d->m_z /= ensembles;
        d->m /= ensembles;
        d->m_stdv = d->m_stdv / ensembles - pow(d->m, 2);
        d->m_stdv = sqrt(max(d->m_stdv, real(0.0)));

        // for debug
//        printf("[naiveAverageMagnetization_kernel] %e, %e, %e, %e, %e\n",
//               d->m_x,
//               d->m_y,
//               d->m_z,
//               d->m,
//               d->m_stdv
//        );
    }
}


void CudaMeasurement::measureAverageMagnetization(std::size_t mstep)
{
    if (mstep % *FortranData::avrg_step != 0)
        return;

    cudaStream_t workStream = CudaParallelizationHelper::def.getWorkStream();

    const uint M = emomM.extent(2);
    const dim3 threads = 1024;
    const dim3 blocks = (M + threads.x - 1) / threads.x;

    // My simple kernel for testing
    naiveAverageMagnetization_kernel<<<blocks, threads, 0, workStream>>>(
            emomM,
            mavg_buff_gpu.data() + mavg_count
    );

    mavg_iter(mavg_count++) = static_cast<uint>(mstep);
    // cudaDeviceSynchronize(); // for printing


    if (mavg_count >= *FortranData::avrg_buff)
    {
        mavg_buff_cpu.copy_sync(mavg_buff_gpu);

        for (uint i = 0; i < mavg_count; ++i)
            measurementWriter.write(MeasurementType::AverageMagnetization,
                                    mavg_iter(i),
                                    (real*)&mavg_buff_cpu[i],
                                    sizeof(AverageMagnetizationData) / sizeof(real));


        mavg_buff_gpu.zeros();
        mavg_buff_cpu.zeros();
        mavg_iter.zeros();
        mavg_count = 0;
    }
}


__global__ void naiveBinderCumulantNoEnergy_kernel(const CudaTensor<real, 3> emomM,
                                                   real temp,
                                                   real mub,
                                                   real k_bolt,
                                                   CudaMeasurement::BinderCumulantData* d)
{
    const uint tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= 1)
        return;

    const uint atoms = emomM.extent(1);
    const uint ensembles = emomM.extent(2);

    for (uint k = 0; k < ensembles; ++k)
    {
        real m[3] = {0};

        for (uint i = 0; i < atoms; ++i)
        {
            m[0] += emomM(0, i, k);
            m[1] += emomM(1, i, k);
            m[2] += emomM(2, i, k);
        }

        const real avrgme = norm(3, m) / atoms;
        const real avrgm2 = pow(avrgme, 2);
        const real avrgm4 = pow(avrgm2, 2);
        d->cumuw += 1;

        d->avrgmcum = (d->avrgmcum * d->cumutotw + avrgme * d->cumuw ) / (d->cumutotw + d->cumuw);
        d->avrgm2cum = (d->avrgm2cum * d->cumutotw + avrgm2 * d->cumuw ) / (d->cumutotw + d->cumuw);
        d->avrgm4cum = (d->avrgm4cum * d->cumutotw + avrgm4 * d->cumuw ) / (d->cumutotw + d->cumuw);

        assert(d->avrgm2cum != 0.0); // this is not checked in Fortran
        d->binderc = 1.0 - (d->avrgm4cum / 3.0 / pow(d->avrgm2cum, 2));

        // For T=0, not the proper susceptibility
//        const real extra = (temp == 0.0)? 1.0 : k_bolt * temp;
//        d->pmsusc = (d->avrgm2cum - pow(d->avrgmcum, 2)) * pow(mub, 2) * atoms / (k_bolt * extra);

        if (temp > 0.0)
            d->pmsusc = (d->avrgm2cum - pow(d->avrgmcum, 2)) * pow(mub, 2) * atoms / pow(k_bolt, 2) / temp;
        else
            d->pmsusc = (d->avrgm2cum - pow(d->avrgmcum, 2)) * pow(mub, 2) * atoms / k_bolt;



        d->cumutotw += d->cumuw;
    }

//    printf("[naiveBinderCumulantNoEnergy_kernel] %e\t%e\t%e\t%e\t%e\t\n",
//           d->avrgmcum,
//           d->avrgm2cum,
//           d->avrgm4cum,
//           d->binderc,
//           d->pmsusc
//    );

}

void CudaMeasurement::measureBinderCumulant(std::size_t mstep)
{
    if (mstep % *FortranData::cumu_step != 0)
        return;


    cudaStream_t workStream = CudaParallelizationHelper::def.getWorkStream();
    if (*FortranData::plotenergy != 0)
    {
        assert(false); // not yet implemented
    }
    else
    {
        naiveBinderCumulantNoEnergy_kernel<<<1, 1, 0, workStream>>>(
                emomM,
                *FortranData::temperature,
                *FortranData::mub,
                *FortranData::k_bolt,
                cumu_buff_gpu.data()
        );
    }

    // cudaDeviceSynchronize(); // for printing without delay

    if (cumu_count++ % *FortranData::cumu_buff == 0) // TODO this is wrong
    {
        cumu_buff_cpu.copy_sync(cumu_buff_gpu);

        measurementWriter.write(MeasurementType::BinderCumulant,
                                    mstep,
                                    (real*)cumu_buff_cpu.data(),
                                    9);


        cumu_buff_gpu.zeros();
        cumu_buff_cpu.zeros();
    }

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

        // printf("[pontryagin_no_kernel] %e\n", *pontryagin_no_out);
    }
}

void CudaMeasurement::measureSkyrmionNumber(std::size_t mstep)
{
    if (mstep % *FortranData::skyno_step != 0)
        return;

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
    grad_mom.zeros();
    grad_moments_kernel<<<blocks, threads, 0, workStream>>>(emom, dxyz_vec, dxyz_atom, dxyz_list, grad_mom);
    pontryagin_no_kernel<<<blocks, threads, 0, workStream>>>(
            emomM,
            grad_mom,
            skyno_buff_gpu.data() + skyno_buff_gpu.extent(0) * skyno_count
    );
    skyno_iter(skyno_count++) = static_cast<uint>(mstep);

    // cudaDeviceSynchronize(); // for print debugging only

    if (skyno_count >= *FortranData::skyno_buff)
    {
        skyno_buff_cpu.copy_sync(skyno_buff_gpu);

        for (uint i = 0; i < skyno_count; ++i)
            measurementWriter.write(MeasurementType::SkyrmionNumber,
                                    skyno_iter(i),
                                    &skyno_buff_cpu(0, i),
                                    skyno_buff_cpu.extent(0));

        skyno_buff_gpu.zeros();
        skyno_buff_cpu.zeros();
        skyno_iter.zeros();
        skyno_count = 0;
    }
}

