// A timer class used for timing tasks.
// Copyright (C) 2018 Anders Logg.

#include <iostream>
#include <iomanip>
#include <ctime>

namespace VirtualCity
{

class Timer
{
public:

    Timer(std::string task) : task(task)
    {
        task[0] = toupper(task[0]);
        std::cout << task + "..." << std::endl;
        clockStart = std::clock();
    }

    ~Timer()
    {
        const clock_t clockStop = std::clock();
        const double elapsedTime = (clockStop - clockStart) / CLOCKS_PER_SEC;
        std::cout << std::fixed << std::setprecision(2)
                  << "Elapsed time (CPU): " << elapsedTime
                  << " (" << task << ")" << std::endl;
    }

private:

    std::string task;
    std::clock_t clockStart;

};

}
