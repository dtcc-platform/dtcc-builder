// Copyright (C) 2020 Anders Logg
// Licensed under the MIT License

#ifndef DTCC_TIMER_H
#define DTCC_TIMER_H

#include <ctime>
#include <map>

namespace DTCC
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
      : Name(std::move(name)), autoStart(autoStart)
  {
    if (autoStart)
      Start();
  }

  // Destructor
  ~Timer()
  {
    if (autoStart)
      Stop();
  }

  // Start clock
  void Start()
  {
    clockStart = std::clock();
  }

  // Stop clock
  void Stop()
  {
    // Record current time
    clockStop = std::clock();

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
    return (clockStop - clockStart) / (double) CLOCKS_PER_SEC;
  }

  // Print elapsed time
  void Print() const
  {
    Info("Elapsed time (CPU): " + str(Time()) + " (" + Name + ")");
  }

  // Print report (summary of all timers)
  static void Report(const std::string& title)
  {
    // Create columns
    std::vector<std::vector<std::string>> cols(4);
    cols[0].push_back(title);
    cols[1].push_back("CPU mean");
    cols[2].push_back("CPU total");
    cols[3].push_back("Count");
    for (const auto &it : timings)
    {
      const std::string name = it.first;
      const double time = it.second.first;
      const size_t count = it.second.second;
      cols[0].push_back(name);
      cols[1].push_back(str(time / static_cast<double>(count)));
      cols[2].push_back(str(time));
      cols[3].push_back(str(count));
    }

    // Compute column sizes
    const size_t tabsize = 2;
    std::vector<size_t> colsize(cols.size());
    std::fill(colsize.begin(), colsize.end(), 0);
    for (size_t i = 0; i < cols[0].size(); i++)
      for (size_t j = 0; j < cols.size(); j++)
        colsize[j] = std::max(colsize[j], cols[j][i].size());
    size_t width = 0;
    for (size_t j = 0; j < cols.size(); j++)
      width += colsize[j] + tabsize;
    width -= tabsize;

    // Generate table
    std::string s = "\n";
    for (size_t i = 0; i < cols[0].size(); i++)
    {
      s += "Timer: ";
      for (size_t j = 0; j < cols.size(); j++)
      {
        s += cols[j][i];
        for (size_t k = 0; k < colsize[j] + tabsize - cols[j][i].size(); k++)
          s += " ";
      }
      s += "\n";
      if (i == 0)
      {
        s += "Timer: ";
        for (size_t k = 0; k < width; k++)
          s += "-";
        s += "\n";
      }
    }

    // Print table
    Info(s);
  }

private:
  // True if timer is started and stopped automatically
  bool autoStart{};

  // Time at start
  std::clock_t clockStart{};

  // Time when stopped
  std::clock_t clockStop{};

  // Array of registered timings: Name --> (Total time, Count)
  static std::map<std::string, std::pair<double, size_t>> timings;
};

// Initialize static member
std::map<std::string, std::pair<double, size_t>> Timer::timings;

} // namespace DTCC

#endif
