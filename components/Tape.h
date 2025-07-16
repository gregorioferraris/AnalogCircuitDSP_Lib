// File: Tape.h
#ifndef TAPE_H
#define TAPE_H

#include "Component.h" // Assumiamo che Component.h definisca la classe base Component
#include <string>
#include <random>    // Per la generazione di rumore
#include <cmath>     // Per funzioni matematiche come atan

/**
 * @class Tape
 * @brief Modella un effetto di nastro magnetico semplificato con guadagno, saturazione, rumore e isteresi.
 *
 * Questa classe simula alcune delle caratteristiche principali di un nastro magnetico
 * analogico, come la saturazione (compressione del segnale), l'aggiunta di rumore
 * e un semplice effetto di isteresi.
 */
class Tape : public Component {
public:
    /**
     * @brief Costruttore per la classe Tape.
     * @param name Nome del componente.
     * @param input_node Nome del nodo di ingresso.
     * @param output_node Nome del nodo di uscita.
     * @param ground_node Nome del nodo di massa.
     * @param saturation_level Livello a cui il segnale inizia a saturare (es. in Weber/metro).
     * @param gain Guadagno lineare complessivo del nastro.
     * @param bias_frequency Frequenza di bias (per linearizzazione e riduzione rumore).
     * @param tape_speed Velocità del nastro (m/s).
     * @param noise_level Livello di rumore additivo (RMS).
     * @param hysteresis_factor Fattore per modellare l'isteresi (semplificato).
     */
    Tape(const std::string& name, const std::string& input_node, const std::string& output_node,
         const std::string& ground_node, double saturation_level = 1.0, double gain = 1.0,
         double bias_frequency = 150e3, double tape_speed = 0.38,
         double noise_level = 0.0, double hysteresis_factor = 0.0);

    /**
     * @brief Processa un singolo campione audio attraverso il modello del nastro.
     * @param input_sample Il campione di ingresso da processare.
     * @return Il campione processato.
     */
    double process(double input_sample) override;

    // Metodi getter per accedere ai parametri (opzionale, ma buona pratica)
    double getSaturationLevel() const { return saturation_level_; }
    double getGain() const { return gain_; }
    double getBiasFrequency() const { return bias_frequency_; }
    double getTapeSpeed() const { return tape_speed_; }
    double getNoiseLevel() const { return noise_level_; }
    double getHysteresisFactor() const { return hysteresis_factor_; }

private:
    double saturation_level_;   ///< Livello di saturazione.
    double gain_;               ///< Guadagno lineare.
    double bias_frequency_;     ///< Frequenza di bias.
    double tape_speed_;         ///< Velocità del nastro.
    double noise_level_;        ///< Livello di rumore.
    double hysteresis_factor_;  ///< Fattore di isteresi.

    double last_output_;        ///< Memorizza l'ultimo campione di uscita per l'isteresi.

    std::default_random_engine random_engine_; ///< Motore per la generazione di numeri casuali.
    std::normal_distribution<double> noise_distribution_; ///< Distribuzione normale per il rumore.
};

#endif // TAPE_H
