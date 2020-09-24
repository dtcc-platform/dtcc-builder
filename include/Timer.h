// A timer class used for timing tasks.
// Copyright (C) 2018 Anders Logg.

#include <ctime>
#include <iomanip>
#include <iostream>

namespace DTCC
{

class Timer
{
public:

  // Create timer for task. By default the clock starts when the timer
  // is constructed and the elapsed time is reported when timer goes
  // out of scope.
  Timer(std::string task, bool autoStart=true) : task(task), autoStart(autoStart)
  {
    if (autoStart)
      Start();
  }

  // Destructor
  ~Timer()
  {
    if (autoStart)
    {
      Stop();
      Print();
    }
  }

  // Start clock
  void Start()
  {
    clockStart = std::clock();
  }

  // Stop clock
  void Stop()
  {
    clockStop = std::clock();
  }

  // Return elapsed time
  double Time() const
  {
    return (clockStop - clockStart) / (double) CLOCKS_PER_SEC;
  }

  // Print elapsed time
  void Print()
  {
    std::cout << std::fixed << std::setprecision(2)
              << "Elapsed time (CPU): " << Time() << " (" << task << ")"
              << std::endl;
  }

private:
  std::string task{};
  bool autoStart{};
  std::clock_t clockStart{};
  std::clock_t clockStop{};
};

} // namespace DTCC
