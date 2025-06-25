#include "cudaMeasurement.cuh"
#include "cudaAverageMagnetization.cuh"
#include "cudaBinderCumulant.cuh"

#include "c_helper.h"
#include "stopwatchPool.hpp"
#include "fortranData.hpp"


CudaMeasurement::CudaMeasurement(const CudaTensor<real, 3>& emomM,
                                 const CudaTensor<real, 3>& emom,
                                 const CudaTensor<real, 2>& mmom,
                                 bool alwaysCopy)
: emomM(emomM)
, emom(emom)
, mmom(mmom)
, stopwatch(GlobalStopwatchPool::get("Cuda measurement"))
, alwaysCopy(alwaysCopy)
{
    if (*FortranData::do_avrg == 'Y' && *FortranData::avrg_step > 0)
    {
        measurables.push_back(std::make_unique<AverageMagnetization>(emomM));
        std::cout << "AverageMagnetization observable added" << std::endl;
    }

    if (*FortranData::do_cumu == 'Y' && *FortranData::cumu_step > 0)
    {
        measurables.push_back(std::make_unique<BinderCumulant>(emomM));
        std::cout << "BinderCumulant observable added" << std::endl;
    }
}


void CudaMeasurement::measure(std::size_t mstep)
{
    // Copy required?
    // TODO: fortran_do_measurements seems to only check flags
//    if ( !(alwaysCopy || fortran_do_measurements(mstep)) )
//        return;

    for (const auto& meas : measurables)
    {
        meas->measure(mstep);
    }
}


void CudaMeasurement::flushMeasurements(std::size_t mstep)
{
 // TODO: perhaps if alwaysCopy is false we could queue up measurements from multiple
 // mstep's and then only sync if buffer is full or if flushMeasurements is called
}

