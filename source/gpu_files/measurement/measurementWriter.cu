#include "measurementWriter.cuh"

#include <iomanip>
#include <iostream>
#include <sstream>


void MeasurementWriter::write(MeasurementType mtype, const size_t* iteration, const void* data, size_t N)
{
    auto [it, inserted] = files.try_emplace(mtype);
    std::ofstream& out = it->second;
    if (inserted)
        initFile(mtype, out);

    for (uint i = 0; i < N; ++i)
    {
        out << std::setw(indent) << iteration[i] << ' ';
        switch (mtype)
        {
            case MeasurementType::AverageMagnetization:
                write(*(static_cast<const AverageMagnetizationData*>(data)+i), out);
                break;

            case MeasurementType::BinderCumulant:
                write(*(static_cast<const BinderCumulantData*>(data)+i), out);
                break;

            case MeasurementType::SkyrmionNumber:
                write(*(static_cast<const SkyrmionNumberData*>(data)+i), out);
                break;
        }
        out << '\n';
    }

}


std::string MeasurementWriter::readSimIDFromFile()
{
    // read the simid from the inpsd.dat file
    const std::string filename = "inpsd.dat";
    std::ifstream inputFile(filename);

    if (!inputFile)
        throw std::runtime_error("Could not open file '" + filename + "'");

    std::string line, keyword, value;
    while (std::getline(inputFile, line))
    {
        std::istringstream iss(line);
        if ((iss >> keyword >> value) && keyword == "simid")
            break;
    }

    if (value.empty())
        throw std::runtime_error("simid not found in " + filename);

    return value;
}



void MeasurementWriter::initFile(MeasurementType mtype, std::ofstream& out)
{
    const std::string filename = this->filename(mtype);

    if (out.is_open())
        out.close();

    // open with truncation (overwrite existing file)
    out.open(filename, std::ios::out | std::ios::trunc);

    if (!out)
        throw std::runtime_error("Could not open file '" + filename + "'");

    out << std::setw(indent);
    for (const auto& headerTag : header(mtype))
        out << headerTag << std::setw(width);
    out << '\n';

    out.close();

    // re-open with appending mode for subsequent writes
    out.open(filename, std::ios::out | std::ios::app);

    if (!out)
        throw std::runtime_error("Could not open file '" + filename + "'");

    out << std::scientific << std::uppercase << std::setprecision(fp_precision);
}


std::string MeasurementWriter::filename(MeasurementType mtype)
{
    std::stringstream ss;
    switch (mtype)
    {
        case MeasurementType::AverageMagnetization:
            ss << "averages";
            break;

        case MeasurementType::BinderCumulant:
            ss << "cumulants";
            break;

        case MeasurementType::SkyrmionNumber:
            ss << "sknumber";
            break;
    }
    ss << "." << simID << ".out";
    return ss.str();
}


std::span<const std::string_view> MeasurementWriter::header(MeasurementType mtype)
{
    using namespace std::string_view_literals;

    static constexpr std::string_view avgMag[] = {
        "#Iter"sv, "<M>_x"sv, "<M>_y"sv, "<M>_z"sv, "<M>"sv, "M_{stdv}"sv
    };

    static constexpr std::string_view binder[] = {
        "#Iter"sv, "<M>"sv, "<M^2>"sv, "<M^4>"sv, "U_{Binder}"sv,
        "\\chi"sv, "C_v(tot)"sv, "<E>"sv, "<E_{exc}>"sv, "<E_{lsf}"sv
    };

    static constexpr std::string_view skyrmion[] = {
        "#Iter"sv, "Skx num"sv, "Skx avg"sv, "Skx std"sv
    };

    switch (mtype)
    {
        case MeasurementType::AverageMagnetization:
            return avgMag;

        case MeasurementType::BinderCumulant:
            return binder;

        case MeasurementType::SkyrmionNumber:
            return skyrmion;
    }
}


void MeasurementWriter::write(const AverageMagnetizationData& data, std::ostream& out)
{
    out << std::setw(width) << data.m_x
        << std::setw(width) << data.m_y
        << std::setw(width) << data.m_z
        << std::setw(width) << data.m
        << std::setw(width) << data.m_stdv;
}


void MeasurementWriter::write(const BinderCumulantData& data, std::ostream& out)
{
    out << std::setw(width) << data.avrgmcum
        << std::setw(width) << data.avrgm2cum
        << std::setw(width) << data.avrgm4cum
        << std::setw(width) << data.binderc
        << std::setw(width) << data.pmsusc
        << std::setw(width) << data.cv
        << std::setw(width) << data.avrgecum
        << std::setw(width) << data.avrgetcum
        << std::setw(width) << data.avrgelcum;
}


void MeasurementWriter::write(const SkyrmionNumberData& data, std::ostream& out)
{
    out << std::setw(width) << data.skyno
        << std::setw(width) << data.skyno_avg
        << std::setw(width) << data.skyno_stdv;
}




