#include "gpuMeasurement.hpp"
#include "c_helper.h"
#include "stopwatchPool.hpp"
#include "gpuParallelizationHelper.hpp"
#include <iostream>
#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
#endif
using ParallelizationHelper = GpuParallelizationHelper;
namespace mm = kernels::measurement;

GpuMeasurement::GpuMeasurement(const GpuTensor<real, 3>& emomM,
                                 const GpuTensor<real, 3>& emom,
                                 const GpuTensor<real, 2>& mmom,
                                 bool alwaysCopy)
: emomM(emomM)
, emom(emom)
, mmom(mmom)
, N(emomM.extent(1))
, M(emomM.extent(2))
, NX(128) // TODO: dont hardcode these values, needs to be imported from Fortran
, NY(128)
, NZ(1)
, NT(1)
, workStream( ParallelizationHelperInstance.getWorkStream())
, stopwatch(GlobalStopwatchPool::get("Gpu measurement"))
, do_avrg(*FortranData::do_avrg == 'Y')
, mavg_kernel_threads(256)
, mavg_kernel_blocks(mm::ceil_div(M, mavg_kernel_threads.x))
, do_cumu(*FortranData::do_cumu == 'Y')
, cumu_kernel_threads(32)
, cumu_kernel_blocks(mm::ceil_div(M, cumu_kernel_threads.x))
, do_skyno([](char c) -> SkyrmionMethod {
                switch (c)
                {
                    case 'Y': return SkyrmionMethod::BruteForce;
                    case 'T': return SkyrmionMethod::Triangulation;
                    default: return SkyrmionMethod::None;
                }
        }(*FortranData::do_skyno))
, nsimp(2 * NX * NY * NZ * NT)
, skyno_kernel_threads(128, 1)
, skyno_kernel_blocks([this](SkyrmionMethod method) -> dim3 {
                switch (method)
                {
                    case SkyrmionMethod::BruteForce: return (mm::ceil_div(N, skyno_kernel_threads.x), M);
                    case SkyrmionMethod::Triangulation: return mm::ceil_div(nsimp, skyno_kernel_threads.x);
                    default: return 0;
                }
        }(do_skyno))
{
    if (do_avrg)
    {
        assert(*FortranData::avrg_step > 0 && *FortranData::avrg_buff > 0);

        mavg_buff_gpu.Allocate(*FortranData::avrg_buff);
        mavg_buff_gpu.zeros();

        mavg_buff_cpu.AllocateHost(*FortranData::avrg_buff);
        mavg_buff_cpu.zeros();

        mavg_partial_buff.Allocate(mavg_kernel_blocks.x);

        mavg_iter.AllocateHost(*FortranData::avrg_buff);
        mavg_iter.zeros();

        std::cout << "AverageMagnetization observable added" << std::endl;
    }

    if (do_cumu)
    {
        assert(*FortranData::cumu_step > 0 && *FortranData::cumu_buff > 0);

        cumu_buff_gpu.Allocate(1);
        cumu_buff_gpu.zeros();

        cumu_buff_cpu.AllocateHost(1);
        cumu_buff_cpu.zeros();

        cumu_partial_buff.Allocate(cumu_kernel_blocks.x);

        std::cout << "BinderCumulant observable added" << std::endl;
    }

    if (do_avrg || do_cumu)
    {
        emomMEnsembleSums.Allocate(3, M);
        emomMEnsembleSums.zeros();
    }

    if (do_skyno == SkyrmionMethod::BruteForce || do_skyno == SkyrmionMethod::Triangulation)
    {
        assert(*FortranData::skyno_step > 0 && *FortranData::skyno_buff > 0);

        skyno_buff_gpu.Allocate(*FortranData::skyno_buff);
        skyno_buff_gpu.zeros();

        skyno_buff_cpu.AllocateHost(*FortranData::skyno_buff);
        skyno_buff_cpu.zeros();

        skyno_iter.AllocateHost(*FortranData::skyno_buff);
        skyno_iter.zeros();

        std::cout << "SkyrmionNumber observable added" << std::endl;
    }

    if (do_skyno == SkyrmionMethod::BruteForce)
    {
        Tensor<real, 3> dxyz_vec_fortran(FortranData::dxyz_vec, 3, 26, N);
        dxyz_vec.Allocate(3, 26, N);
        dxyz_vec.copy_sync(dxyz_vec_fortran);

        Tensor<int, 2> dxyz_atom_fortran(FortranData::dxyz_atom, 26, N);
        dxyz_atom.Allocate(26, N);
        dxyz_atom.copy_sync(dxyz_atom_fortran);

        Tensor<int, 1> dxyz_list_fortran(FortranData::dxyz_list, N);
        dxyz_list.Allocate(N);
        dxyz_list.copy_sync(dxyz_list_fortran);

        grad_mom.Allocate(3, 3, N, M);
        grad_mom.zeros();

        skyno_partial_buff.Allocate(skyno_kernel_threads.x * skyno_kernel_threads.y);
    }

    if (do_skyno == SkyrmionMethod::Triangulation)
    {
        simp.Allocate(3, nsimp);
        skyno_partial_buff.Allocate(skyno_kernel_blocks.x);

        const uint pairs = nsimp / 2;
        const uint threads = 256;
        const uint blocks = mm::ceil_div(pairs, threads);

        mm::delaunay_tri_tri<<<blocks, threads, 0, workStream>>>(NX, NY, NZ, NT, simp);
    }

    stopwatch.add("constructor");
}


GpuMeasurement::~GpuMeasurement()
{
    if (do_avrg)
    {
        mavg_buff_gpu.Free();
        mavg_buff_cpu.FreeHost();
        mavg_partial_buff.Free();
        mavg_iter.FreeHost();
    }

    if (do_cumu)
    {
        cumu_buff_gpu.Free();
        cumu_buff_cpu.FreeHost();
        cumu_partial_buff.Free();
    }

    if (do_avrg || do_cumu)
    {
        emomMEnsembleSums.Free();
    }

    if (do_skyno == SkyrmionMethod::BruteForce || do_skyno == SkyrmionMethod::Triangulation)
    {
        skyno_buff_gpu.Free();
        skyno_buff_cpu.FreeHost();
        skyno_iter.FreeHost();
        skyno_partial_buff.Free();
    }

    if (do_skyno == SkyrmionMethod::BruteForce)
    {
        dxyz_vec.Free();
        dxyz_atom.Free();
        dxyz_list.Free();
        grad_mom.Free();
    }

    if (do_skyno == SkyrmionMethod::Triangulation)
    {
        simp.Free();
    }
}


void GpuMeasurement::measure(std::size_t mstep)
{
    --mstep; // this is because the simulation loop begins at 1 because of Fortran indexing

    const bool avrg = timeToMeasure(MeasurementType::AverageMagnetization, mstep);
    const bool cumu = timeToMeasure(MeasurementType::BinderCumulant, mstep);

    if (avrg || cumu)
    {
        calculateEmomMSum();
        stopwatch.add("sum reduction of emomM for shared use");
    }

    if (avrg)
    {
        measureAverageMagnetization(mstep);
        stopwatch.add("average magnetization");
    }

    if (cumu)
    {
        measureBinderCumulant(mstep);
        stopwatch.add("binder cumulant");
    }

    if (timeToMeasure(MeasurementType::SkyrmionNumber, mstep))
    {
        measureSkyrmionNumber(mstep);
        stopwatch.add("skyrmion number");
    }
}


void GpuMeasurement::flushMeasurements(std::size_t mstep)
{
    if (do_avrg)
        saveToFile(MeasurementType::AverageMagnetization);

    if (do_cumu)
        saveToFile(MeasurementType::BinderCumulant);

    if (do_skyno != SkyrmionMethod::None)
        saveToFile(MeasurementType::SkyrmionNumber);
}


void GpuMeasurement::measureAverageMagnetization(std::size_t mstep)
{
    const size_t smem = mm::nwarps(mavg_kernel_threads) * sizeof(real);

    mm::averageMagnetization_partial<<<mavg_kernel_blocks, mavg_kernel_threads, smem, workStream>>>(
            emomMEnsembleSums, N, M, mavg_partial_buff.data()
    );

    mm::averageMagnetization_finalize<<<1, mavg_kernel_threads, smem, workStream>>>(
            mavg_partial_buff.data(), mavg_kernel_blocks.x, M, mavg_buff_gpu.data()[mavg_count]
    );

    mavg_iter(mavg_count++) = mstep;

    if (mavg_count >= *FortranData::avrg_buff)
    {
        saveToFile(MeasurementType::AverageMagnetization);
    }
}


void GpuMeasurement::measureBinderCumulant(std::size_t mstep)
{
    if (*FortranData::plotenergy == 0)
    {
        const size_t smem = mm::nwarps(cumu_kernel_threads) * sizeof(real);

        mm::binderCumulantNoEnergy_partial<<<cumu_kernel_blocks, cumu_kernel_threads, smem, workStream>>>(
                emomMEnsembleSums, N, M, cumu_partial_buff.data()
        );

        mm::binderCumulantNoEnergy_finalize<<<1, cumu_kernel_threads, smem, workStream>>>(
                cumu_partial_buff.data(),
                cumu_kernel_blocks.x,
                N,
                M,
                *FortranData::temperature,
                *FortranData::mub,
                *FortranData::k_bolt,
                *cumu_buff_gpu.data()
        );
    }
    else
    {
        throw std::invalid_argument("Not yet implemented.");
    }


    if ((cumu_count++ % *FortranData::cumu_buff) == 0)
    {
        saveToFile(MeasurementType::BinderCumulant);
    }

}


void GpuMeasurement::measureSkyrmionNumber(std::size_t mstep)
{
    if (do_skyno == SkyrmionMethod::BruteForce)
    {
        mm::grad_moments<<<skyno_kernel_blocks, skyno_kernel_threads, 0, workStream>>>(
                emomM, dxyz_vec, dxyz_atom, dxyz_list, grad_mom
        );


        size_t smem = mm::nwarps(skyno_kernel_threads) * sizeof(real);
        mm::pontryagin_no_partial<<<skyno_kernel_blocks, skyno_kernel_threads, smem, workStream>>>(
                emomM, grad_mom, skyno_partial_buff.data()
        );

        smem = skyno_kernel_threads.x * sizeof(real);
        mm::pontryagin_no_finalize<<<1, skyno_kernel_threads, smem, workStream>>>(
                skyno_partial_buff.data(), skyno_kernel_blocks.x, M, skyno_count + 1, skyno_buff_gpu.data()[skyno_count]
        );
    }

    else if (do_skyno == SkyrmionMethod::Triangulation)
    {
        size_t smem = mm::nwarps(skyno_kernel_threads) * sizeof(real);
        mm::pontryagin_tri_partial<<<skyno_kernel_blocks, skyno_kernel_threads, smem, workStream>>>(
                emom, simp, skyno_partial_buff.data()
        );

        smem = skyno_kernel_threads.x * sizeof(real);
        mm::pontryagin_tri_finalize<<<1, skyno_kernel_threads, smem, workStream>>>(
                skyno_partial_buff.data(), skyno_kernel_blocks.x, M, skyno_count + 1, skyno_buff_gpu.data()[skyno_count]
        );
    }

    skyno_iter(skyno_count++) = mstep;

    if (skyno_count >= *FortranData::skyno_buff)
    {
        saveToFile(MeasurementType::SkyrmionNumber);
    }
}


void GpuMeasurement::calculateEmomMSum()
{
    emomMEnsembleSums.zeros();
    const dim3 threads(256, 1);
    const dim3 blocks(mm::ceil_div(N, threads.x), 3 * M);
    const size_t smem = mm::nwarps(threads) * sizeof(real);
    mm::sumOverAtoms<<<blocks, threads, smem, workStream>>>(emomM, emomMEnsembleSums);
}


void GpuMeasurement::saveToFile(MeasurementType mtype)
{
    switch (mtype)
    {
        case MeasurementType::AverageMagnetization:
            mavg_buff_cpu.copy_sync(mavg_buff_gpu);

            measurementWriter.write(
                    mavg_iter.data(),
                    mavg_buff_cpu.data(),
                    mavg_count
            );

            mavg_buff_gpu.zeros();
            mavg_buff_cpu.zeros();
            mavg_iter.zeros();
            mavg_count = 0;
        break;

        case MeasurementType::BinderCumulant:
            cumu_buff_cpu.copy_sync(cumu_buff_gpu);

            measurementWriter.write(
                    &cumu_count, // TODO: in fortran the equiv to cumu_count is printed, but should it not be mstep?
                    cumu_buff_cpu.data(),
                    1
            );
        break;

        case MeasurementType::SkyrmionNumber:
            skyno_buff_cpu.copy_sync(skyno_buff_gpu);

            measurementWriter.write(
                    skyno_iter.data(),
                    skyno_buff_cpu.data(),
                    skyno_count
            );

            skyno_buff_gpu.zeros();
            skyno_buff_cpu.zeros();
            skyno_iter.zeros();
            skyno_count = 0;
        break;

        default:
            throw std::invalid_argument("Not yet implemented.");
    }
}


bool GpuMeasurement::timeToMeasure(MeasurementType mtype, size_t mstep) const
{
    switch (mtype)
    {
        case MeasurementType::AverageMagnetization:
            return do_avrg && ((mstep % *FortranData::avrg_step) == 0);

        case MeasurementType::BinderCumulant:
            return do_cumu && ((mstep % *FortranData::cumu_step) == 0);

        case MeasurementType::SkyrmionNumber:
            return do_skyno != SkyrmionMethod::None && ((mstep % *FortranData::skyno_step) == 0);

        default:
            throw std::invalid_argument("Not yet implemented.");
    }
}

