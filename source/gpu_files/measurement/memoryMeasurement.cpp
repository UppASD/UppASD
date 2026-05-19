#include "memoryMeasurement.h"


void TensorMemoryTracker::printResults()
{
    const auto get_unit = [](size_t x) -> std::pair<double, std::string_view> {
        if (x >= GB) return {static_cast<double>(GB), "GB"};
        else if (x >= MB) return {static_cast<double>(MB), "MB"};
        else if (x >= KB) return {static_cast<double>(KB), "KB"};
        else return {1.0, "B"};
    };

    auto& t = get();
    const std::scoped_lock lock(t.host_mutex, t.device_mutex);

    const auto max_host = *std::max_element(t.host_memory.cbegin(), t.host_memory.cend());
    const auto max_device = *std::max_element(t.device_memory.cbegin(), t.device_memory.cend());

    const auto [h_divisor, h_unit] = get_unit(max_host);
    const auto [d_divisor, d_unit] = get_unit(max_device);

    std::cout << "\n------------------TENSOR MEMORY CONSUMPTION REPORT----------------------\n";

    std::cout << std::setfill(' ') << std::left
              << std::setw(colWidth) << "CATEGORY" << std::setw(colWidth) << "PEAK CONSUMPTION" << '\n'
              << std::setw(colWidth) << "Max CPU"
              << std::setw(colWidth) << (static_cast<double>(max_host)   / h_divisor) << ' ' << h_unit << '\n'
              << std::setw(colWidth) << "Max GPU"
              << std::setw(colWidth) << (static_cast<double>(max_device) / d_divisor)<< ' ' << d_unit << '\n';

    std::cout << "-----------------END CPU-GPU DATA MOVEMENT REPORT-------------------\n";
}


void TensorMemoryTracker::saveToFile()
{
    auto& t = get();
    const std::scoped_lock lock(t.host_mutex, t.device_mutex);

    std::ofstream host_file("memory_over_time_host.csv", std::ios::out | std::ios::trunc);
    if (!host_file) throw std::runtime_error("Could not open file 'memory_over_time_host.csv'");

    std::ofstream device_file("memory_over_time_device.csv", std::ios::out | std::ios::trunc);
    if (!device_file) throw std::runtime_error("Could not open file 'memory_over_time_device.csv'");

    host_file << "bytes\n";

    for (const auto& mem : t.host_memory)
        host_file << mem << '\n';

    device_file << "bytes\n";
    for (const auto& mem : t.device_memory)
        device_file << mem << '\n';
}


void TensorMemoryTracker::add_host(int64_t bytes)
{
    const std::lock_guard<std::mutex> lock(get().host_mutex);
    auto& mem = get().host_memory;
    const int64_t prev = mem.back();
    mem.push_back(prev + bytes);
}


void TensorMemoryTracker::add_device(int64_t bytes)
{
    const std::lock_guard<std::mutex> lock(get().device_mutex);
    auto& mem = get().device_memory;
    const int64_t prev = mem.back();
    mem.push_back(prev + bytes);
}


void TensorMemoryTracker::remove_host(int64_t bytes)
{
    const std::lock_guard<std::mutex> lock(get().host_mutex);
    auto& mem = get().host_memory;
    const int64_t prev = mem.back();
    mem.push_back(prev - bytes);
}


void TensorMemoryTracker::remove_device(int64_t bytes)
{
    const std::lock_guard<std::mutex> lock(get().device_mutex);
    auto& mem = get().device_memory;
    const int64_t prev = mem.back();
    mem.push_back(prev - bytes);
}


void TensorDataMovementTracker::printResults()
{
    const auto& t = get();

    const uint64_t h2d       = t.h2d.load(std::memory_order_relaxed);
    const uint64_t d2h       = t.d2h.load(std::memory_order_relaxed);
    const uint64_t h2d_async = t.h2d_async.load(std::memory_order_relaxed);
    const uint64_t d2h_async = t.d2h_async.load(std::memory_order_relaxed);

    const auto any_bigger = [&](uint64_t thr) -> bool {
        return h2d >= thr || d2h >= thr || h2d_async >= thr || d2h_async >= thr;
    };

    double divisor = 1.0;
    std::string_view unit = "B";
    if (any_bigger(GB))      { divisor = static_cast<double>(GB); unit = "GB"; }
    else if (any_bigger(MB)) { divisor = static_cast<double>(MB); unit = "MB"; }
    else if (any_bigger(KB)) { divisor = static_cast<double>(KB); unit = "kB"; }

    print(divisor, unit, h2d, h2d_async, d2h, d2h_async);
}


// all add methods are thread safe
void TensorDataMovementTracker::add_h2d(size_t bytes)       { get().h2d.fetch_add(bytes, std::memory_order_relaxed); }
void TensorDataMovementTracker::add_d2h(size_t bytes)       { get().d2h.fetch_add(bytes, std::memory_order_relaxed); }
void TensorDataMovementTracker::add_h2d_async(size_t bytes) { get().h2d_async.fetch_add(bytes, std::memory_order_relaxed); }
void TensorDataMovementTracker::add_d2h_async(size_t bytes) { get().d2h_async.fetch_add(bytes, std::memory_order_relaxed); }


void TensorDataMovementTracker::print(double divisor, std::string_view unit,
                  uint64_t h2d, uint64_t h2d_async,
                  uint64_t d2h, uint64_t d2h_async)
{

    std::cout << "\n------------------CPU-GPU DATA MOVEMENT REPORT----------------------\n";

    std::cout << std::setfill(' ') << std::left
              << std::setw(colWidth) << "CATEGORY"         << std::setw(colWidth) << unit << '\n'
              << std::setw(colWidth) << "CPU -> GPU"       << std::setw(colWidth) << (static_cast<double>(h2d)       / divisor) << '\n'
              << std::setw(colWidth) << "CPU -> GPU Async" << std::setw(colWidth) << (static_cast<double>(h2d_async) / divisor) << '\n'
              << std::setw(colWidth) << "GPU -> CPU"       << std::setw(colWidth) << (static_cast<double>(d2h)       / divisor) << '\n'
              << std::setw(colWidth) << "GPU -> CPU Async" << std::setw(colWidth) << (static_cast<double>(d2h_async) / divisor) << '\n';

    std::cout << "-----------------END CPU-GPU DATA MOVEMENT REPORT-------------------\n";
}
