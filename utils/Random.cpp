// utils/Random.cpp
#include "Random.h"
#include <random> // Per std::mt19937, std::uniform_real_distribution
#include <chrono> // Per std::chrono::system_clock::now()
#include <iostream> // Per debug

// Inizializzazione del generatore di numeri casuali
// Usiamo un seed basato sul tempo per una maggiore casualità ad ogni esecuzione
static std::mt19937 rng(std::chrono::system_clock::now().time_since_epoch().count());

// Funzione per generare un numero casuale in un intervallo [min, max]
double getRandomDouble(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

// Funzione per applicare una tolleranza percentuale a un valore nominale
// Restituisce un valore casuale all'interno di +/- tolerance_percent dal valore nominale.
double applyTolerance(double nominal_value, double tolerance_percent) {
    if (tolerance_percent < 0.0) {
        // La tolleranza non può essere negativa, la trattiamo come 0.
        // std::cerr << "Warning: tolerance_percent cannot be negative. Using 0% tolerance." << std::endl;
        return nominal_value;
    }
    double deviation = nominal_value * (tolerance_percent / 100.0);
    return getRandomDouble(nominal_value - deviation, nominal_value + deviation);
}
