#include "cudaMeasurement.cuh"
#include "cudaAverageMagnetization.cuh"
#include "cudaBinderCumulant.cuh"
#include "cudaSkyrmionNumber.cuh"

#include "c_helper.h"
#include "stopwatchPool.hpp"
#include "fortranData.hpp"


CudaMeasurement::CudaMeasurement(const CudaTensor<real, 3>& emomM,
                                 const CudaTensor<real, 3>& emom,
                                 const CudaTensor<real, 2>& mmom,
                                 bool alwaysCopy)
: stopwatch(GlobalStopwatchPool::get("Cuda measurement"))
, alwaysCopy(alwaysCopy)
{
    if (*FortranData::do_avrg == 'Y')
    {
        assert(*FortranData::avrg_step > 0 && *FortranData::avrg_buff > 0);
        measurables.push_back(std::make_unique<AverageMagnetization>(emomM));
        std::cout << "AverageMagnetization observable added" << std::endl;
    }

    if (*FortranData::do_cumu == 'Y')
    {
        assert(*FortranData::cumu_step > 0 && *FortranData::cumu_buff > 0);
        measurables.push_back(std::make_unique<BinderCumulant>(emomM));
        std::cout << "BinderCumulant observable added" << std::endl;
    }

    if (*FortranData::do_skyno == 'Y')
    {
        assert(*FortranData::skyno_step > 0 && *FortranData::skyno_buff > 0);
        measurables.push_back(std::make_unique<SkyrmionNumber>(emomM, emom));
        std::cout << "SkyrmionNumber observable added" << std::endl;
    }
}


void CudaMeasurement::measure(std::size_t mstep)
{
    // Copy required?
    // TODO: fortran_do_measurements seems to only check flags
//    if ( !(alwaysCopy || fortran_do_measurements(mstep)) )
//        return;
    --mstep;
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

