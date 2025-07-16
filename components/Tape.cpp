// File: Tape.cpp
#include "Tape.h"
#include <iostream> // Per output di debug, se necessario

/**
 * @brief Implementazione del costruttore della classe Tape.
 * @param name Nome del componente.
 * @param input_node Nome del nodo di ingresso.
 * @param output_node Nome del nodo di uscita.
 * @param ground_node Nome del nodo di massa.
 * @param saturation_level Livello a cui il segnale inizia a saturare.
 * @param gain Guadagno lineare complessivo del nastro.
 * @param bias_frequency Frequenza di bias.
 * @param tape_speed Velocità del nastro.
 * @param noise_level Livello di rumore additivo.
 * @param hysteresis_factor Fattore per modellare l'isteresi.
 */
Tape::Tape(const std::string& name, const std::string& input_node, const std::string& output_node,
           const std::string& ground_node, double saturation_level, double gain,
           double bias_frequency, double tape_speed, double noise_level,
           double hysteresis_factor)
    : Component(name, input_node, output_node, ground_node), // Chiama il costruttore della classe base
      saturation_level_(saturation_level),
      gain_(gain),
      bias_frequency_(bias_frequency),
      tape_speed_(tape_speed),
      noise_level_(noise_level),
      hysteresis_factor_(hysteresis_factor),
      last_output_(0.0), // Inizializza l'ultimo output a zero
      random_engine_(std::random_device{}()), // Inizializza il motore casuale
      noise_distribution_(0.0, 1.0) // Distribuzione normale con media 0 e deviazione standard 1
{
    // Potresti aggiungere qui logica di inizializzazione specifica
    // ad esempio, validazione dei parametri.
    // std::cout << "Tape '" << name_ << "' creato." << std::endl;
}

/**
 * @brief Processa un singolo campione audio attraverso il modello del nastro.
 *
 * L'elaborazione include:
 * 1. Applicazione del guadagno lineare.
 * 2. Saturazione tramite funzione atan (approssimazione morbida).
 * 3. Aggiunta di rumore gaussiano.
 * 4. Applicazione di una semplice isteresi (feedback del campione precedente).
 * @param input_sample Il campione di ingresso da processare.
 * @return Il campione processato.
 */
double Tape::process(double input_sample) {
    // 1. Applica il guadagno lineare
    double processed_sample = input_sample * gain_;

    // 2. Saturazione (usando atan per una curva di saturazione morbida)
    // Normalizza il segnale rispetto al livello di saturazione prima di applicare atan
    if (saturation_level_ > 0) {
        processed_sample = saturation_level_ * std::atan(processed_sample / saturation_level_);
    }

    // 3. Aggiungi rumore
    if (noise_level_ > 0) {
        processed_sample += noise_distribution_(random_engine_) * noise_level_;
    }

    // 4. Isteresi semplificata: una combinazione lineare dell'input processato
    // e dell'output precedente, pesata dal fattore di isteresi.
    // Questo crea una leggera "memoria" del segnale.
    processed_sample = (1.0 - hysteresis_factor_) * processed_sample +
                       hysteresis_factor_ * last_output_;

    // Aggiorna l'ultimo output per il prossimo ciclo
    last_output_ = processed_sample;

    return processed_sample;
}
