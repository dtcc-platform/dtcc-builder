// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_TIMER_H
#define DTCC_TIMER_H

#include <chrono>
#include <map>
#include <utility>

// #include "JSON.h"
#include "Table.h"

namespace DTCC_BUILDER
{

class Timer
{
public:
  /// Name of timer
  std::string Name{};

  /// Create timer. By default the clock starts when the timer
  /// is constructed and the elapsed time is reported when timer
  /// goes out of scope.
  explicit Timer(std::string name, bool autoStart = true)
      : Name(std::move(name)), autoStart(autoStart), running(false)
  {
    if (autoStart)
      Start();
  }

  // Destructor
  ~Timer()
  {
    if (autoStart and running)
      Stop();
  }

  // Start clock
  void Start()
  {
    t0 = std::chrono::high_resolution_clock::now();
    running = true;
  }

  // Stop clock
  void Stop()
  {
    // Check if timer is running
    if (!running)
    {
      warning("Timer stopped but it's not running");
      return;
    }

    // Record current time
    t1 = std::chrono::high_resolution_clock::now();
    running = false;

    // Register timing
    auto it = timings.find(Name);
    if (it == timings.end())
    {
      timings[Name] = std::make_pair(Time(), 1);
    }
    else
    {
      it->second.first += Time();
      it->second.second += 1;
    }
  }

  // Return elapsed time
  double Time() const
  {
    std::chrono::duration<double, std::milli> dt = t1 - t0;
    return dt.count() / 1000.0;
  }

  // Print elapsed time
  void Print() const
  {
    info("Elapsed time (CPU): " + str(Time()) + " (" + Name + ")");
  }

  // Print report (summary of all timers)
  static void Report(const std::string &title,
                     const std::string &fileName = "Timings.json")
  {
    // Build table
    Table table(title);
    table.Rows.push_back({"Task", "CPU mean", "CPU total", "Count"});
    for (const auto &it : timings)
    {
      const std::string task = it.first;
      const double total = it.second.first;
      const size_t count = it.second.second;
      const double mean = total / static_cast<double>(count);
      std::vector<std::string> row;
      row.push_back(task);
      row.push_back(str(mean));
      row.push_back(str(total));
      row.push_back(str(count));
      table.Rows.push_back(row);
    }

    // Print table
    info(table);

    // Write JSON
    // JSON::Write(timings, fileName, 4);
  }

private:
  // True if timer is started and stopped automatically
  bool autoStart{};

  // True if timer is running
  bool running{};

  // Time at start
  std::chrono::high_resolution_clock::time_point t0{};

  // Time when stopped
  std::chrono::high_resolution_clock::time_point t1{};

  // Array of registered timings: Name --> (Total time, Count)
  static std::map<std::string, std::pair<double, size_t>> timings;
};

// Initialize static member
std::map<std::string, std::pair<double, size_t>> Timer::timings;

} // namespace DTCC_BUILDER

#endif
