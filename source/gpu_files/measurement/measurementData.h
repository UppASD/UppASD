#pragma once

#include <cstdint>
#include "real_type.h"
#include "fortranData.hpp"

enum class MeasurementType : uint8_t
{
    AverageMagnetization, BinderCumulant, SkyrmionNumber
};

enum class SkyrmionMethod : uint8_t
{
    None,        // no measurement
    BruteForce,  // gradient + Pontryagin density method
    Triangulation
};


//struct MeasurementConfig
//{
//    // Average magnetization
//    bool do_avrg{};
//    size_t avrg_step{};
//    size_t avrg_buff{};
//
//    // Binder cumulant
//    bool do_cumu{};
//    size_t cumu_step{};
//    size_t cumu_buff{};
//
//    // Skyrmion
//    SkyrmionMethod skyno_method{SkyrmionMethod::None};
//    size_t skyno_step{};
//    size_t skyno_buff{};
//
//    // build this once from FortranData in a small helper
//    static MeasurementConfig makeConfigFromFortran();
//};


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
    real cv{};            // Specific heat
    real avrgecum{};      // Cumulated average of E
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

