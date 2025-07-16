// utils/Logger.cpp
#include "Logger.h"
#include <iostream>
#include <iomanip> // Per std::fixed e std::setprecision
#include <chrono> // Per timestamp
#include <ctime> // Per std::localtime, std::mktime
#include <sstream> // Per std::stringstream

// Inizializzazione della variabile statica del livello di log
LogLevel Logger::currentLogLevel = LogLevel::INFO;

// Funzione per impostare il livello di log
void Logger::setLogLevel(LogLevel level) {
    currentLogLevel = level;
}

// Funzione per ottenere il timestamp corrente formattato
std::string Logger::getTimestamp() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    // Converti in tm struct per la formattazione
    std::tm bt = *std::localtime(&in_time_t);

    // Aggiungi i millisecondi
    auto duration = now.time_since_epoch();
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(duration) % 1000;

    std::stringstream ss;
    ss << std::put_time(&bt, "%Y-%m-%d %H:%M:%S");
    ss << '.' << std::setfill('0') << std::setw(3) << milliseconds.count();
    return ss.str();
}

// Funzione di log generica
void Logger::log(LogLevel level, const std::string& message) {
    if (level >= currentLogLevel) {
        std::string levelStr;
        switch (level) {
            case LogLevel::DEBUG: levelStr = "DEBUG"; break;
            case LogLevel::INFO:  levelStr = "INFO "; break; // Spazio per allineamento
            case LogLevel::WARN:  levelStr = "WARN "; break; // Spazio per allineamento
            case LogLevel::ERROR: levelStr = "ERROR"; break;
            case LogLevel::FATAL: levelStr = "FATAL"; break;
        }
        std::cout << "[" << getTimestamp() << "] [" << levelStr << "] " << message << std::endl;
    }
}
