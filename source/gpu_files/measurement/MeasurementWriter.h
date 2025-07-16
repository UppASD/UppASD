#pragma once

#include <map>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include "tensor.cuh"
#include "real_type.h"


enum class MeasurementType
{
    AverageMagnetization, BinderCumulant, SkyrmionNumber
};


class MeasurementWriter
{
public:
    MeasurementWriter();
    void write(MeasurementType mtype, size_t iteration, const real* data, uint N);

private:
    std::map<MeasurementType, std::ofstream> outputFiles;

    std::string simID;

    std::string outputFilename(MeasurementType mtype);
    static std::vector<std::string> fileHeader(MeasurementType mtype);
};
