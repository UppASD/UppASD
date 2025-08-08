#pragma once

#include "real_type.h"


enum class MeasurementType
{
    AverageMagnetization, BinderCumulant, SkyrmionNumber
};

struct AverageMagnetizationData
{
    real m_x, m_y, m_z, m, m_stdv;
};


struct BinderCumulantData
{
    real avrgmcum;      // Cumulated average of m
    real avrgm2cum;     // Cumulated average of m^2
    real avrgm4cum;     // Cumulated average of m^4
    real binderc;       // Binder cumulant
    real pmsusc;        // Susceptibility
    real cv;            // Specific heat
    real avrgecum;      // Cumulated average of E
    // real avrge2cum;     // Cumulated average of E^2
    real avrgetcum;     // Cumulated average of E_xc
    real avrgelcum;     // Cumulated average of E_LSF

    real cumuw;         // Weight for current sample to cumulant
    real cumutotw;      // Sum of all cumulant weights
    uint Navrgcum;      // Counter for number of cumulated averages
};


struct SkyrmionNumberData
{
    real skyno, skyno_avg, skyno_stdv;
};

