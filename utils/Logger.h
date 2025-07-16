// utils/Logger.h
#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip> // For std::put_time
#include <ctime>   // For std::localtime

// Enum for log levels
enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR,
    FATAL
};

// Logger class (Singleton for easy global access)
class Logger {
public:
    // Get the Logger instance
    static Logger& getInstance();

    // Initialize the logger (called once at application/plugin startup)
    void init(LogLevel minLevel = LogLevel::INFO, const std::string& filename = "");

    // Methods to log messages at different levels
    void log(LogLevel level, const std::string& message,
             const std::string& file = "", int line = 0, const std::string& func = "");

    // Helper functions for different log levels
    void debug(const std::string& message, const std::string& file = "", int line = 0, const std::string& func = "");
    void info(const std::string& message, const std::string& file = "", int line = 0, const std::string& func = "");
    void warning(const std::string& message, const std::string& file = "", int line = 0, const std::string& func = "");
    void error(const std::string& message, const std::string& file = "", int line = 0, const std::string& func = "");
    void fatal(const std::string& message, const std::string& file = "", int line = 0, const std::string& func = "");

private:
    Logger(); // Private constructor for Singleton pattern
    ~Logger(); // Destructor

    // Delete copy constructor and assignment operator
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    std::ofstream logFile;
    LogLevel minLogLevel;
    bool initialized;

    std::string getTimestamp();
    std::string levelToString(LogLevel level);
};

// Macros for convenient logging (automatically include file, line, function)
#define LOG_DEBUG(msg) Logger::getInstance().debug(msg, __FILE__, __LINE__, __func__)
#define LOG_INFO(msg) Logger::getInstance().info(msg, __FILE__, __LINE__, __func__)
#define LOG_WARNING(msg) Logger::getInstance().warning(msg, __FILE__, __LINE__, __func__)
#define LOG_ERROR(msg) Logger::getInstance().error(msg, __FILE__, __LINE__, __func__)
#define LOG_FATAL(msg) Logger::getInstance().fatal(msg, __FILE__, __LINE__, __func__)

#endif // LOGGER_H
