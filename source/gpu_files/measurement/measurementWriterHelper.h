#pragma once

#include <concepts>
#include <ostream>
#include <array>
#include <string_view>
#include <iomanip>

#include "measurementData.h"

// Per data type metadata
template<class T>
struct MeasurementTraits;


template<class T>
concept MeasurementTypeLike = requires(const T& t, std::ostream& os, int w)
{
    { MeasurementTraits<T>::filebase } -> std::convertible_to<std::string_view>;
    MeasurementTraits<T>::columns;
    { MeasurementTraits<T>::print(t, os, w) } -> std::same_as<void>;
};


template<class... Ts>
void write_cols(std::ostream& out, int width, const Ts&... xs)
{
    ((out << std::setw(width) << xs), ...);
}


template<> struct MeasurementTraits<AverageMagnetizationData>
{
    static constexpr std::string_view filebase = "averages";

    static constexpr std::array<std::string_view, 6> columns = {
        "#Iter", "<M>_x", "<M>_y", "<M>_z", "<M>", "M_{stdv}"
    };

    static void print(const AverageMagnetizationData& a, std::ostream& out, int width)
    {
        write_cols(out, width, a.m_x, a.m_y, a.m_z, a.m, a.m_stdv);
    }
};


template<> struct MeasurementTraits<BinderCumulantData>
{
    static constexpr std::string_view filebase = "cumulants";

    static constexpr std::array<std::string_view, 10> columns = {
        "#Iter", "<M>", "<M^2>", "<M^4>", "U_{Binder}",
        "\\chi", "C_v(tot)", "<E>", "<E_{exc}>", "<E_{lsf}>"
    };

    static void print(const BinderCumulantData& b, std::ostream& out, int width)
    {
        write_cols(out, width,
                   b.avrgmcum, b.avrgm2cum, b.avrgm4cum, b.binderc,
                   b.pmsusc, b.cv, b.avrgecum, b.avrgetcum, b.avrgelcum);
    }
};


template<> struct MeasurementTraits<SkyrmionNumberData>
{
    static constexpr std::string_view filebase = "sknumber";

    static constexpr std::array<std::string_view, 4> columns = {
        "#Iter", "Skx num", "Skx avg", "Skx std"
    };

    static void print(const SkyrmionNumberData& s, std::ostream& out, int width)
    {
        write_cols(out, width, s.skyno, s.skyno_avg, s.skyno_stdv);
    }
};

