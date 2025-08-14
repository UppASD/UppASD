#include "cudaMeasurement.cuh"

#include <cuda_runtime.h>

#include "c_helper.h"
#include "stopwatchPool.hpp"

#include "cudaParallelizationHelper.hpp"
#include "kernels.cuh"

#include <iostream>


CudaMeasurement::CudaMeasurement(const CudaTensor<real, 3>& emomM,
                                 const CudaTensor<real, 3>& emom,
                                 const CudaTensor<real, 2>& mmom,
                                 bool alwaysCopy)
: emomM(emomM)
, emom(emom)
, mmom(mmom)
, N(emomM.extent(1))
, M(emomM.extent(2))
, workStream(CudaParallelizationHelper::def.getWorkStream())
, stopwatch(GlobalStopwatchPool::get("Cuda measurement"))
, do_avrg(*FortranData::do_avrg == 'Y')
, do_cumu(*FortranData::do_cumu == 'Y')
, do_skyno([](char c) -> SkyrmionMethod {
                switch (c)
                {
                    case 'Y': return SkyrmionMethod::BruteForce;
                    case 'T': return SkyrmionMethod::Triangulation;
                    default: return SkyrmionMethod::None;
                }
        }(*FortranData::do_skyno))
, nsimp()
{
    if (do_avrg)
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

    if (do_cumu)
    {
        assert(*FortranData::cumu_step > 0 && *FortranData::cumu_buff > 0);

        cumu_buff_gpu.Allocate(1);
        cumu_buff_gpu.zeros();

        cumu_buff_cpu.AllocateHost(1);
        cumu_buff_cpu.zeros();

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
    }

    if (do_skyno == SkyrmionMethod::Triangulation)
    {
        simp.Allocate(3, nsimp);
    }

    stopwatch.add("constructor");
}


CudaMeasurement::~CudaMeasurement()
{
    if (do_avrg)
    {
        mavg_buff_gpu.Free();
        mavg_buff_cpu.FreeHost();
        mavg_iter.FreeHost();
    }

    if (do_cumu)
    {
        cumu_buff_gpu.Free();
        cumu_buff_cpu.FreeHost();
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


void CudaMeasurement::measure(std::size_t mstep)
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


void CudaMeasurement::flushMeasurements(std::size_t mstep)
{
    if (do_avrg)
        saveToFile(MeasurementType::AverageMagnetization);

    if (do_cumu)
        saveToFile(MeasurementType::BinderCumulant);

    if (do_skyno != SkyrmionMethod::None)
        saveToFile(MeasurementType::SkyrmionNumber);
}


void CudaMeasurement::measureAverageMagnetization(std::size_t mstep)
{
    const uint threads = std::min(256u, M);
    const uint blocks = (M + threads - 1) / threads;

    kernels::averageMagnetization<<<blocks, threads, 0, workStream>>>(
        emomMEnsembleSums,
        N,
        M,
        mavg_buff_gpu.data()[mavg_count]
    );

    mavg_iter(mavg_count++) = mstep;

    if (mavg_count >= *FortranData::avrg_buff)
    {
        saveToFile(MeasurementType::AverageMagnetization);
    }
}



void CudaMeasurement::measureBinderCumulant(std::size_t mstep)
{
    if (*FortranData::plotenergy == 0)
    {
        kernels::binderCumulantNoEnergy<<<1, 1, 0, workStream>>>(
            emomMEnsembleSums,
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
        assert(false); // not yet implemented
    }


    if ((cumu_count++ % *FortranData::cumu_buff) == 0)
    {
        saveToFile(MeasurementType::BinderCumulant);
    }

}


void CudaMeasurement::measureSkyrmionNumber(std::size_t mstep)
{
    constexpr dim3 threads = 256;
    const dim3 blocks = {
        (N + threads.x - 1) / threads.x,
        (M + threads.y - 1) / threads.y
    };

    // TODO are we sure this should be emom, and not emomM?
    kernels::grad_moments<<<blocks, threads, 0, workStream>>>(emom, dxyz_vec, dxyz_atom, dxyz_list, grad_mom);
    kernels::pontryagin_no<<<blocks, threads, 0, workStream>>>(
        emomM,
        grad_mom,
        skyno_count + 1,
        skyno_buff_gpu.data()[skyno_count]
    );

    skyno_iter(skyno_count++) = mstep;

    if (skyno_count >= *FortranData::skyno_buff)
    {
        saveToFile(MeasurementType::SkyrmionNumber);
    }
}


void CudaMeasurement::calculateEmomMSum()
{
    constexpr dim3 threads = 256;

    dim3 blocks = {3, M};
    kernels::sumOverAtoms<<<blocks, threads, 0, workStream>>>(emomM, emomMEnsembleSums);
}


void CudaMeasurement::saveToFile(MeasurementType mtype)
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


bool CudaMeasurement::timeToMeasure(MeasurementType mtype, size_t mstep) const
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

