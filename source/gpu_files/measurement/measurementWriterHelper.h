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
concept MeasurementTypeLike = requires(const T& t, std::ostream& os, int w, bool flag)
{
    { MeasurementTraits<T>::filebase } -> std::convertible_to<std::string_view>;
    { MeasurementTraits<T>::write_header(os, w, flag) } -> std::same_as<void>;
    { MeasurementTraits<T>::print(t, os, w, flag) } -> std::same_as<void>;
};


template<class... Ts>
void write_cols(std::ostream& out, int width, const Ts&... xs)
{
    ((out << std::setw(width) << xs), ...);
}


template<> struct MeasurementTraits<AverageMagnetizationData>
{
    static constexpr std::string_view filebase = "averages";

    static void write_header(std::ostream& out, int w, bool)
    {
        write_cols(out, w, "#Iter","<M>_x","<M>_y","<M>_z","<M>","M_{stdv}");
    }

    static void print(const AverageMagnetizationData& a, std::ostream& out, int w, bool)
    {
        write_cols(out, w, a.m_x, a.m_y, a.m_z, a.m, a.m_stdv);
    }
};


template<> struct MeasurementTraits<BinderCumulantData>
{
    static constexpr std::string_view filebase = "cumulants";

    static void write_header(std::ostream& out, int w, bool)
    {
        write_cols(out, w,
                   "#Iter","<M>","<M^2>","<M^4>","U_{Binder}",
                   "\\chi","C_v(tot)","<E>","<E_{exc}>","<E_{lsf}>");
    }

    static void print(const BinderCumulantData& b, std::ostream& out, int w, bool)
    {
        write_cols(out, w,
                   b.avrgmcum, b.avrgm2cum, b.avrgm4cum, b.binderc,
                   b.pmsusc, b.cv, b.avrgecum, b.avrgetcum, b.avrgelcum);
    }
};


template<> struct MeasurementTraits<SkyrmionNumberData>
{
    static constexpr std::string_view filebase = "sknumber";

    static void write_header(std::ostream& out, int w, bool)
    {
        write_cols(out, w, "#Iter","Skx num","Skx avg","Skx std");
    }

    static void print(const SkyrmionNumberData& s, std::ostream& out, int w, bool)
    {
        write_cols(out, w, s.skyno, s.skyno_avg, s.skyno_stdv);
    }
};



template<> struct MeasurementTraits<EnergyData>
{
    static constexpr std::string_view filebase = "totenergy";

    static void write_header(std::ostream& out, int w, bool do_jtensor)
    {
        if (do_jtensor)
        {
            write_cols(out, w,
                "#Iter","Tot","Heis-Tens","Ani","PD","BiqDM","BQ",
                "Dip","Zeeman","LSF","Chir","Ring","SA");
        }
        else
        {
            write_cols(out, w,
                "#Iter","Tot","Exc","Ani","DM","PD","BiqDM","BQ",
                "Dip","Zeeman","LSF","Chir","Ring","SA");
        }
    }

    static void print(const EnergyData& d, std::ostream& out, int w, bool do_jtensor)
    {
        if (do_jtensor)
        {
            write_cols(out, w,
                d.total,
                d.pair,
                d.anisotropy,
                d.PD,
                d.BiqDM,
                d.BQ,
                d.Dip,
                d.Zeeman,
                d.LSF,
                d.Chir,
                d.Ring,
                d.SA);
        }
        else
        {
            write_cols(out, w,
                d.total,
                d.exchange,
                d.anisotropy,
                d.DM,
                d.PD,
                d.BiqDM,
                d.BQ,
                d.Dip,
                d.Zeeman,
                d.LSF,
                d.Chir,
                d.Ring,
                d.SA);
        }
    }
};


template<> struct MeasurementTraits<EnergyStdData>
{
    static constexpr std::string_view filebase = "stdenergy";

    static void write_header(std::ostream& out, int w, bool do_jtensor)
    {
        if (do_jtensor)
        {
            write_cols(out, w,
                "#Iter","Tot","Heis-Tens","Ani","PD","BiqDM","BQ",
                "Dip","Zeeman","LSF","Chir","Ring","SA");
        }
        else
        {
            write_cols(out, w,
                "#Iter","Tot","Exc","Ani","DM","PD","BiqDM","BQ",
                "Dip","Zeeman","LSF","Chir","Ring","SA");
        }
    }

    static void print(const EnergyStdData& v, std::ostream& out, int w, bool do_jtensor)
    {
        const auto& d = v.ene_ref;

        if (do_jtensor)
        {
            write_cols(out, w,
                d.std_total,
                d.std_pair,
                d.std_anisotropy,
                d.std_PD,
                d.std_BiqDM,
                d.std_BQ,
                d.std_Dip,
                d.std_Zeeman,
                d.std_LSF,
                d.std_Chir,
                d.std_Ring,
                d.std_SA);
        }
        else
        {
            write_cols(out, w,
                d.std_total,
                d.std_exchange,
                d.std_anisotropy,
                d.std_DM,
                d.std_PD,
                d.std_BiqDM,
                d.std_BQ,
                d.std_Dip,
                d.std_Zeeman,
                d.std_LSF,
                d.std_Chir,
                d.std_Ring,
                d.std_SA);
        }
    }
};



template<>
struct MeasurementTraits<AutocorrelationData>
{
    static constexpr std::string_view filebase = "autocorr";

    static void write_header(std::ostream&, int, bool) {}

    static void print(const AutocorrelationData& a, std::ostream& out, int width, bool)
    {
        for (size_t j = 0; j < a.size; ++j)
            out << std::setw(width) << a.values[j];
    }
};