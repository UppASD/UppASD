#pragma once

#include <cstddef>
#include <fstream>
#include <ios>
#include <iomanip>
#include <string>
#include <typeindex>
#include <unordered_map>
#include <vector>
#include <ostream>

#include "measurementData.h"
#include "measurementWriterHelper.h"


class MeasurementWriter
{
public:
    explicit MeasurementWriter(int fp_precision = 8, int padding = 6);

    template<class Iter_t, MeasurementTypeLike Data_t>
    void write(const Iter_t* iteration, const Data_t* data, size_t N);

private:
    template<MeasurementTypeLike T>
    std::ofstream& getFile();

    template<MeasurementTypeLike T>
    void initFile(std::ofstream& out) const;

    template<MeasurementTypeLike T>
    std::string filename() const;

    static std::string readSimIDFromFile();

private:
    static constexpr int fp_printed_symbols = 6;
    const int fp_precision;
    const int colWidth;

    const std::string simID = readSimIDFromFile();
    std::unordered_map<std::type_index, std::ofstream> files;
};


template<class Iter_t, MeasurementTypeLike Data_t>
void MeasurementWriter::write(const Iter_t* iteration, const Data_t* data, size_t N)
{
    std::ofstream& out = getFile<Data_t>();
    for (size_t i = 0; i < N; ++i)
    {
        out << std::setw(colWidth) << iteration[i];
        MeasurementTraits<Data_t>::print(data[i], out, colWidth);
        out << '\n';
    }
}


template<MeasurementTypeLike T>
std::ofstream& MeasurementWriter::getFile()
{
    auto [it, inserted] = files.try_emplace(std::type_index(typeid(T)));
    if (inserted)
        initFile<T>(it->second);
    return it->second;
}


template<MeasurementTypeLike T>
void MeasurementWriter::initFile(std::ofstream& out) const
{
    const std::string fname = filename<T>();
    out.open(fname, std::ios::out | std::ios::trunc);
    if (!out)
        throw std::runtime_error("Could not open file '" + fname + "'");

    out << std::setfill(' ') << std::right;

    for (const auto& h : MeasurementTraits<T>::columns)
        out << std::setw(colWidth) << h;
    out << '\n';

    out << std::scientific << std::uppercase << std::setprecision(fp_precision);
}


template<MeasurementTypeLike T>
std::string MeasurementWriter::filename() const
{
    std::string s;
    s.reserve(MeasurementTraits<T>::filebase.size() + 1 + simID.size() + 4);
    s.append(MeasurementTraits<T>::filebase);
    s.push_back('.');
    s.append(simID);
    s.append(".out");
    return s;
}

