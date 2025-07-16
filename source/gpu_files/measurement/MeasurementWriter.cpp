#include "MeasurementWriter.h"

#include <iomanip>
#include <iostream>


MeasurementWriter::MeasurementWriter()
{
    // read the simid from the inpsd.dat file
    const std::string filename = "inpsd.dat";
    std::ifstream inputFile(filename);
    if (!inputFile.is_open())
        throw std::runtime_error("Could not open input file: " + filename);

    std::string line;
    while (std::getline(inputFile, line))
    {
        std::istringstream iss(line);
        std::string keyword, value;
        if ((iss >> keyword >> value) && keyword == "simid")
        {
            simID = value;
            break;
        }
    }

    if (simID.empty())
        throw std::runtime_error("simid not found in " + filename);
}


void MeasurementWriter::write(MeasurementType mtype, size_t iteration, const real* data, uint N)
{
    const auto indent = std::setw(8);
    const auto colWidth = std::setw(16);

    const std::string filename = outputFilename(mtype);

    if (iteration == 0)
    {
        std::cout << "[MeasurementWriter::write] iteration == 0, filename = " << filename << std::endl;
        // can be checked with "contains" if C++20 is allowed
        auto it = outputFiles.find(mtype);
        if (it != outputFiles.end() && it->second.is_open())
            it->second.close();

        // Open with truncation (overwrite)
        std::ofstream& out = outputFiles[mtype];
        out.open(filename, std::ios::out | std::ios::trunc);

        if (!out)
            throw std::runtime_error("Failed to open file '" + filename + "' for writing.");

        out << indent;
        for (const auto& headerTag : fileHeader(mtype))
            out << headerTag << colWidth;
        out << '\n';

        out.close();
        if (!out)
            throw std::runtime_error("Failed to close file '" + filename + "' for writing.");

        out.open(filename, std::ios::out | std::ios::app);
        if (!out)
            throw std::runtime_error("Failed to open file '" + filename + "' for writing.");

        out << std::scientific << std::uppercase << std::setprecision(8);
    }

    std::ofstream& out = outputFiles[mtype];
    out << indent << iteration;
    for (uint i = 0; i < N; ++i)
        out << colWidth << data[i];
    out << "\n";

    if (!out)
        throw std::runtime_error("Failed to write to file '" + filename + "'.");
}


std::string MeasurementWriter::outputFilename(MeasurementType mtype)
{
    std::stringstream ss;
    ss << "test.";
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

        default:
            throw std::invalid_argument("MeasurementType enum value "
                + std::to_string(static_cast<int>(mtype)) + " not yet implemented.");
    }
    ss << '.' << simID << ".out";
    return ss.str();
}


std::vector<std::string> MeasurementWriter::fileHeader(MeasurementType mtype)
{
    switch (mtype)
    {
        case MeasurementType::AverageMagnetization:
            return {"#Iter", "<M>_x", "<M>_y", "<M>_z", "<M>", "M_{stdv}"};

        case MeasurementType::BinderCumulant:
            return {"#Iter", "<M>", "<M^2>", "<M^4>", "U_{Binder}", "\\chi", "C_v(tot)", "<E>", "<E_{exc}>", "<E_{lsf}"};

        case MeasurementType::SkyrmionNumber:
            return {"#Iter", "Skx num", "Skx avg", "Skx std"};

        default:
            throw std::invalid_argument("MeasurementType enum value "
                                        + std::to_string(static_cast<int>(mtype)) + " not yet implemented.");
    }
}





