#pragma once

#include <cstdint>
#include "real_type.h"
#include "fortranData.hpp"

enum class MeasurementType : uint8_t
{
    AverageMagnetization, BinderCumulant, SkyrmionNumber, Energy, Autocorrelation
};

enum class SkyrmionMethod : uint8_t
{
    None,        // no measurement
    BruteForce,  // gradient + Pontryagin density method
    Triangulation
};

inline SkyrmionMethod skyrmionMethodFromFlag(char c)
{
    switch (c)
    {
        case 'Y': return SkyrmionMethod::BruteForce;
        case 'T': return SkyrmionMethod::Triangulation;
        default: return SkyrmionMethod::None;
    }
}

struct AverageMagnetizationData
{
    real m_x{};
    real m_y{};
    real m_z{};
    real m{};
    real m_stdv{};
};

struct BinderCumulantData
{
    real avrgmcum{};      // Cumulated average of m
    real avrgm2cum{};     // Cumulated average of m^2
    real avrgm4cum{};     // Cumulated average of m^4
    real binderc{};       // Binder cumulant
    real pmsusc{};        // Susceptibility
    real avrgecum{};      // Cumulated average of E
    real cv{};            // Specific heat

    // real avrge2cum;     // Cumulated average of E^2
    real avrgetcum{};     // Cumulated average of E_xc
    real avrgelcum{};     // Cumulated average of E_LSF

    real cumuw{};         // Weight for current sample to cumulant
    real cumutotw{};      // Sum of all cumulant weights
    uint Navrgcum{};      // Counter for number of cumulated averages
};

struct SkyrmionNumberData
{
    real skyno{};
    real skyno_avg{};
    real skyno_stdv{};
};

struct EnergyData
{
    real total;
    real exchange;
    real anisotropy;
    real DM;
    real Zeeman;

    real std_total;
    real std_exchange;
    real std_anisotropy;
    real std_DM;
    real std_Zeeman;


    // fields below not yet implemented
    real PD, BiqDM, BQ, Dip, LSF, Chir, Ring, SA;
    real std_PD, std_BiqDM, std_BQ, std_Dip, std_LSF, std_Chir, std_Ring, std_SA;

};
struct AutocorrelationData
{
    const real* values{};
    size_t size{};
};
