// JSON I/O
// Anders Logg 2019

#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <string>

namespace VirtualCity
{

class CommandLine
{
public:

    static std::string GetOption(std::string option, int argc, char* argv[])
    {
        auto end = argv + argc;
        auto it = find(argv, end, option);
        if (it != end && ++it != end)
            return *it;
        return 0;
    }

};

}

#endif
