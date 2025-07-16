// File: Oscillator.h
#ifndef OSCILLATOR_H
#define OSCILLATOR_H

#include <cmath>   // For std::sin, std::fmod
#include <string>  // For std::string (if needed for debugging/naming)

/**
 * @class Oscillator
 * @brief Genera vari tipi di forme d'onda (sinusoidale, quadra, triangolare, dente di sega).
 *
 * Questa classe è un generatore di segnali autonomo, utile per la sintesi
 * audio o come sorgente di segnale in simulazioni.
 */
class Oscillator {
public:
    /**
     * @brief Enum per i diversi tipi di forme d'onda.
     */
    enum WaveformType {
        SINE,       ///< Onda sinusoidale
        SQUARE,     ///< Onda quadra
        TRIANGLE,   ///< Onda triangolare
        SAWTOOTH    ///< Onda a dente di sega
    };

    /**
     * @brief Costruttore per la classe Oscillator.
     * @param sample_rate La frequenza di campionamento in Hz (es. 44100.0).
     * @param frequency La frequenza dell'onda in Hz (es. 440.0 per La4).
     * @param amplitude L'ampiezza di picco dell'onda (es. 1.0 per il massimo).
     * @param type Il tipo di forma d'onda da generare (default SINE).
     */
    Oscillator(double sample_rate, double frequency = 440.0,
               double amplitude = 1.0, WaveformType type = SINE);

    /**
     * @brief Genera il prossimo campione della forma d'onda.
     * @return Il valore del campione generato.
     */
    double process();

    // Metodi setter per i parametri
    void setFrequency(double freq);
    void setAmplitude(double amp);
    void setWaveformType(WaveformType type);
    void setSampleRate(double rate);
    void resetPhase(); // Resetta la fase corrente a zero

    // Metodi getter
    double getFrequency() const { return frequency_; }
    double getAmplitude() const { return amplitude_; }
    WaveformType getWaveformType() const { return type_; }
    double getSampleRate() const { return sample_rate_; }

private:
    WaveformType type_;         ///< Tipo di forma d'onda corrente.
    double sample_rate_;        ///< Frequenza di campionamento in Hz.
    double frequency_;          ///< Frequenza dell'onda in Hz.
    double amplitude_;          ///< Ampiezza di picco dell'onda.
    double current_phase_;      ///< Fase corrente dell'onda in radianti (da 0 a 2*PI).
    double phase_increment_;    ///< Incremento di fase per ogni campione.

    /**
     * @brief Calcola l'incremento di fase basato su frequenza e sample rate.
     */
    void calculatePhaseIncrement();
};

#endif // OSCILLATOR_H
