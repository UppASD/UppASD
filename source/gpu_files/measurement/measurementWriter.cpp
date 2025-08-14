#include "measurementWriter.h"
#include <fstream>
#include <sstream>
#include <stdexcept>


MeasurementWriter::MeasurementWriter(int fp_precision, int padding)
: fp_precision(fp_precision)
, colWidth(fp_precision + padding + fp_printed_symbols)
{
    if (padding <= 0)
        throw std::invalid_argument("Padding must be bigger than zero.");
}


std::string MeasurementWriter::readSimIDFromFile()
{
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
