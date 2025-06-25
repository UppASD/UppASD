#pragma once

#include "measurable.hpp"

struct BinderCumulantData
{
    real avrgmcum;      // Cumulated average of m
    real avrgm2cum;     // Cumulated average of m^2
    real avrgm4cum;     // Cumulated average of m^4
    real binderc;       // Binder cumulant
    real pmsusc;        // Susceptibility
    real cv;            // Specific heat
    real avrgecum;      // Cumulated average of E
    real avrge2cum;     // Cumulated average of E^2
    real avrgetcum;     // Cumulated average of E_xc
    real avrgelcum;     // Cumulated average of E_LSF
    real cumuw;         // Weight for current sample to cumulant
    real cumutotw;      // Sum of all cumulant weights
    uint Navrgcum;      // Counter for number of cumulated averages

    // these are used for calc_and_print_afm_cumulant:
    // real avrgl2cum;     // Cumulated average of l^2
    // real avrgl4cum;     // Cumulated average of l^4
};

class BinderCumulant : public Measurable
{
public:
    BinderCumulant(const CudaTensor<real, 3>& emomM,
                   dim3 cudaThreads = 1024);
    ~BinderCumulant() override;
    // need ene%energy, ene%ene_xc and ene_lsf for plotenergy
    void measure(std::size_t mstep) override;
    void flushMeasurements(std::size_t mstep) override;
    // void copyToFortran() override;

private:
    const CudaTensor<real, 3>& emomM; // not really needed if AverageMagnetization is injected
    // const AverageMagnetization& mavg; // could be injected, but then we get coupling and cannot run cumu without avrg
    const uint cumu_buff_size;
    const uint N;
    const uint M;
    const uint cumu_step;
    const real temp;
    const bool plotenergy;
    const dim3 threads;
    const dim3 blocks;

    CudaTensor<BinderCumulantData, 1> cumu_buff;
    uint cumu_buff_count = 0;


    
    // TODO: add Tensor for Fortran data pointer
};
