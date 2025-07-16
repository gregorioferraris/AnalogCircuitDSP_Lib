// utils/Random.h
#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <chrono>

// Funzione per generare un numero casuale in un intervallo [min, max]
inline double getRandomDouble(double min, double max) {
    // Usa un generatore di numeri casuali basato sul tempo per una maggiore casualità
    static std::mt19937_64 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

// Funzione per applicare una tolleranza percentuale a un valore nominale
inline double applyTolerance(double nominal_value, double tolerance_percent) {
    if (tolerance_percent == 0.0) {
        return nominal_value;
    }
    double min_val = nominal_value * (1.0 - tolerance_percent / 100.0);
    double max_val = nominal_value * (1.0 + tolerance_percent / 100.0);
    return getRandomDouble(min_val, max_val);
}

#endif // RANDOM_H
