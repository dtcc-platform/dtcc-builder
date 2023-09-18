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
  /// name of timer
  std::string name{};

  /// Create timer. By default the clock starts when the timer
  /// is constructed and the elapsed time is reported when timer
  /// goes out of scope.
  explicit Timer(std::string name, bool auto_start = true)
      : name(std::move(name)), auto_start(auto_start), running(false)
  {
    if (auto_start)
      start();
  }

  // Destructor
  ~Timer()
  {
    if (auto_start and running)
      stop();
  }

  // start clock
  void start()
  {
    t0 = std::chrono::high_resolution_clock::now();
    running = true;
  }

  // stop clock
  void stop()
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
    auto it = timings.find(name);
    if (it == timings.end())
    {
      timings[name] = std::make_pair(time(), 1);
    }
    else
    {
      it->second.first += time();
      it->second.second += 1;
    }
  }

  // Return elapsed time
  double time() const
  {
    std::chrono::duration<double, std::milli> dt = t1 - t0;
    return dt.count() / 1000.0;
  }

  // print elapsed time
  void print() const
  {
    info("Elapsed time (CPU): " + str(time()) + " (" + name + ")");
  }

  // print report (summary of all timers)
  static void report(const std::string &title,
                     const std::string &file_name = "Timings.json")
  {
    // build table
    Table table(title);
    table.rows.push_back({"Task", "CPU mean", "CPU total", "Count"});
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
      table.rows.push_back(row);
    }

    // print table
    info(table);

    // Write JSON
    // JSON::Write(timings, file_name, 4);
  }

private:
  // True if timer is started and stopped automatically
  bool auto_start{};

  // True if timer is running
  bool running{};

  // time at start
  std::chrono::high_resolution_clock::time_point t0{};

  // time when stopped
  std::chrono::high_resolution_clock::time_point t1{};

  // Array of registered timings: name --> (Total time, Count)
  static std::map<std::string, std::pair<double, size_t>> timings;
};

// Initialize static member
std::map<std::string, std::pair<double, size_t>> Timer::timings;

} // namespace DTCC_BUILDER

#endif
