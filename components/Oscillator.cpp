// File: Oscillator.cpp
#include "Oscillator.h"
#include <cmath> // Per M_PI (o definisci PI se non disponibile)

// Definisci PI se M_PI non è disponibile in cmath su alcuni compilatori
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Costruttore per la classe Oscillator.
 * @param sample_rate La frequenza di campionamento in Hz.
 * @param frequency La frequenza dell'onda in Hz.
 * @param amplitude L'ampiezza di picco dell'onda.
 * @param type Il tipo di forma d'onda da generare.
 */
Oscillator::Oscillator(double sample_rate, double frequency,
                       double amplitude, WaveformType type)
    : type_(type),
      sample_rate_(sample_rate),
      frequency_(frequency),
      amplitude_(amplitude),
      current_phase_(0.0),
      phase_increment_(0.0)
{
    // Assicurati che i parametri siano validi
    if (sample_rate_ <= 0) sample_rate_ = 44100.0;
    if (frequency_ < 0) frequency_ = 0.0;
    if (amplitude_ < 0) amplitude_ = 0.0; // L'ampiezza non può essere negativa

    calculatePhaseIncrement();
}

/**
 * @brief Calcola l'incremento di fase basato su frequenza e sample rate.
 * Questo metodo viene chiamato ogni volta che la frequenza o il sample rate cambiano.
 */
void Oscillator::calculatePhaseIncrement() {
    phase_increment_ = (2.0 * M_PI * frequency_) / sample_rate_;
}

/**
 * @brief Genera il prossimo campione della forma d'onda.
 *
 * Questo metodo calcola il valore del campione in base al tipo di forma d'onda
 * e aggiorna la fase corrente per il prossimo campione.
 *
 * @return Il valore del campione generato.
 */
double Oscillator::process() {
    double sample = 0.0;

    switch (type_) {
        case SINE:
            sample = amplitude_ * std::sin(current_phase_);
            break;
        case SQUARE:
            // L'onda quadra è positiva per la prima metà del ciclo, negativa per la seconda.
            // Si basa sul segno della funzione seno per semplicità.
            sample = (std::sin(current_phase_) >= 0) ? amplitude_ : -amplitude_;
            break;
        case TRIANGLE: {
            // Normalizza la fase da 0 a 1 per facilitare i calcoli
            double normalized_phase = current_phase_ / (2.0 * M_PI);
            // Mappa la fase normalizzata al range -1.0 a 1.0 per l'onda triangolare
            if (normalized_phase < 0.25) { // Primo quarto
                sample = amplitude_ * (4.0 * normalized_phase);
            } else if (normalized_phase < 0.75) { // Secondo e terzo quarto
                sample = amplitude_ * (1.0 - 4.0 * (normalized_phase - 0.25));
            } else { // Ultimo quarto
                sample = amplitude_ * (-1.0 + 4.0 * (normalized_phase - 0.75));
            }
            break;
        }
        case SAWTOOTH: {
            // Normalizza la fase da 0 a 1
            double normalized_phase = current_phase_ / (2.0 * M_PI);
            // Mappa la fase normalizzata al range -1.0 a 1.0 per l'onda a dente di sega
            sample = amplitude_ * (2.0 * normalized_phase - 1.0);
            break;
        }
        default:
            // Non dovrebbe succedere, ma come fallback
            sample = 0.0;
            break;
    }

    // Aggiorna la fase per il prossimo campione
    current_phase_ += phase_increment_;
    // Mantiene la fase all'interno dell'intervallo [0, 2*PI)
    current_phase_ = std::fmod(current_phase_, 2.0 * M_PI);
    // Assicurati che la fase sia sempre positiva (fmod può restituire negativi)
    if (current_phase_ < 0) {
        current_phase_ += 2.0 * M_PI;
    }

    return sample;
}

// Implementazioni dei metodi setter
void Oscillator::setFrequency(double freq) {
    if (freq >= 0) {
        frequency_ = freq;
        calculatePhaseIncrement(); // Ricalcola l'incremento di fase
    }
}

void Oscillator::setAmplitude(double amp) {
    if (amp >= 0) {
        amplitude_ = amp;
    }
}

void Oscillator::setWaveformType(WaveformType type) {
    type_ = type;
}

void Oscillator::setSampleRate(double rate) {
    if (rate > 0) {
        sample_rate_ = rate;
        calculatePhaseIncrement(); // Ricalcola l'incremento di fase
    }
}

void Oscillator::resetPhase() {
    current_phase_ = 0.0;
}
