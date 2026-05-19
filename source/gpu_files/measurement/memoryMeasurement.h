#pragma once

#include <iostream>
#include <iomanip>
#include <atomic>
#include <mutex>
#include <cstdint>
#include <cstddef>
#include <string_view>
#include <vector>
#include <algorithm>
#include <fstream>

class TensorMemoryTracker
{
public:
    TensorMemoryTracker(const TensorMemoryTracker&) = delete;
    TensorMemoryTracker& operator=(const TensorMemoryTracker&) = delete;
    static void printResults();
    static void saveToFile();
    static void add_host(int64_t bytes);
    static void add_device(int64_t bytes);
    static void remove_host(int64_t bytes);
    static void remove_device(int64_t bytes);

private:
    TensorMemoryTracker() = default;
    static TensorMemoryTracker& get() { static TensorMemoryTracker m{}; return m; }

private:
    static constexpr int colWidth = 24;
    static constexpr int64_t KB = 1e3;
    static constexpr int64_t MB = 1e6;
    static constexpr int64_t GB = 1e9;
    std::mutex host_mutex;
    std::mutex device_mutex;
    std::vector<int64_t> host_memory{0};
    std::vector<int64_t> device_memory{0};
};


class TensorDataMovementTracker
{
public:
    TensorDataMovementTracker(const TensorDataMovementTracker&) = delete;
    TensorDataMovementTracker& operator=(const TensorDataMovementTracker&) = delete;
    static void printResults();
    static void add_h2d(size_t bytes);
    static void add_d2h(size_t bytes);
    static void add_h2d_async(size_t bytes);
    static void add_d2h_async(size_t bytes);

private:
    TensorDataMovementTracker() = default;
    static TensorDataMovementTracker& get() { static TensorDataMovementTracker m{}; return m; }
    static void print(double divisor, std::string_view unit,
                      uint64_t h2d, uint64_t h2d_async,
                      uint64_t d2h, uint64_t d2h_async);

private:
    static constexpr int colWidth = 24;
    static constexpr uint64_t KB = 1e3;
    static constexpr uint64_t MB = 1e6;
    static constexpr uint64_t GB = 1e9;
    std::atomic<uint64_t> h2d{0};
    std::atomic<uint64_t> d2h{0};
    std::atomic<uint64_t> h2d_async{0};
    std::atomic<uint64_t> d2h_async{0};
};