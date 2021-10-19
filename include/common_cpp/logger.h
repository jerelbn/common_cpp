#pragma once

#include <fstream>

namespace common
{

class Logger
{
public:
    Logger() {}

    Logger(const std::string& filename)
    {
        open(filename);
    }

    void open(const std::string& filename)
    {
        file_.open(filename);
    }

    ~Logger()
    {
        file_.close();
    }

    template <typename... T>
    void log(const T&... data)
    {
        (file_.write(reinterpret_cast<const char*>(&data), sizeof(T)), ...);
    }

private:
    std::ofstream file_;
};


} // namespace common
