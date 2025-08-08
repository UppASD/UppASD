#pragma once

#include <unordered_map>
#include <fstream>
#include <string>
#include <string_view>
#include <span>
#include "measurementData.h"


class MeasurementWriter
{
public:
    void write(MeasurementType mtype, const size_t* iteration, const void* data, size_t N);

private:
    void initFile(MeasurementType mtype, std::ofstream& out);
    std::string readSimIDFromFile();
    std::span<const std::string_view> header(MeasurementType mtype);
    std::string filename(MeasurementType mtype);

    void write(const AverageMagnetizationData& data, std::ostream& out);
    void write(const BinderCumulantData& data, std::ostream& out);
    void write(const SkyrmionNumberData& data, std::ostream& out);

private:
    constexpr static int indent = 8;
    constexpr static int width = 16;
    constexpr static int fp_precision = 8;
    const std::string simID = readSimIDFromFile();

    std::unordered_map<MeasurementType, std::ofstream> files;
};

