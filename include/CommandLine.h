// JSON I/O
// Anders Logg 2019

#ifndef VC_COMMAND_LINE_H
#define VC_COMMAND_LINE_H

#include <string>

namespace VirtualCity
{

class CommandLine
{
public:

    static bool HasOption(std::string option, int argc, char* argv[])
    {
        for (size_t i = 1; i < argc; i++)
            if (option == argv[i])
                return true;
        return false;
    }

    static bool EndsWith(const std::string string, const std::string& ending)
    {
        if (ending.size() > string.size())
            return false;
        return std::equal(ending.rbegin(), ending.rend(), string.rbegin());
    }

    static std::string GetOption(std::string option, int argc, char* argv[])
    {
        for (size_t i = 1; i < argc; i++)
            if (option == argv[i])
                return argv[i + 1];
        return "";
    }

    static int GetIntOption(std::string option, int argc, char* argv[])
    {
        return atoi(GetOption(option, argc, argv).c_str());
    }

};

}

#endif
