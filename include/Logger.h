#include <iostream>

class Logger
{
public:
    static void log (const std::string message)
    {
        std::cout << "[LogInfo] " << message << std::endl;
    }
    static void logError (const std::string message)
    {
        std::cerr << "[LogError]" << message << std::endl;
    }
    static void logWarning (const std::string message)
    {
        std::cout << "[LogWarning] " << message << std::endl;
    }
};