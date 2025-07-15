#include "cudaEnergy.cuh"
#include "fortranData.hpp"

Energy::Energy()
: avrg_step(*FortranData::avrg_step)
, buffer_size(*FortranData::eavrg_buff)
{
    const uint M = 1;
    eavg_buff.Allocate(3, buffer_size, M);
}


Energy::~Energy()
{
    eavg_buff.Free();
}


void Energy::measure(std::size_t mstep)
{
    if ((mstep-1) % avrg_step != 0)
        return;

    // exc=exc+update_ene(emomM(1:3,ii,kk),beff_xc,0.5_dblprec)
}


void Energy::flushMeasurements(std::size_t mstep)
{

}
