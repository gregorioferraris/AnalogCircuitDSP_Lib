// File: Pentode.cpp
#include "Pentode.h"
#include <algorithm> // Per std::max

/**
 * @brief Implementazione del costruttore della classe Pentode.
 * @param name Nome del componente.
 * @param input_node Nome del nodo di ingresso.
 * @param output_node Nome del nodo di uscita.
 * @param ground_node Nome del nodo di massa.
 * @param gain Fattore di guadagno lineare.
 * @param saturation_level Livello di ingresso a cui inizia la saturazione.
 * @param cutoff_threshold Livello di ingresso sotto il quale il segnale viene tagliato.
 * @param output_dc_offset Offset DC aggiunto all'uscita.
 * @param knee_factor Fattore che influenza la "durezza" della curva di saturazione.
 */
Pentode::Pentode(const std::string& name, const std::string& input_node, const std::string& output_node,
                 const std::string& ground_node, double gain, double saturation_level,
                 double cutoff_threshold, double output_dc_offset, double knee_factor)
    : Component(name, input_node, output_node, ground_node), // Chiama il costruttore della classe base
      gain_(gain),
      saturation_level_(saturation_level),
      cutoff_threshold_(cutoff_threshold),
      output_dc_offset_(output_dc_offset),
      knee_factor_(knee_factor)
{
    // Assicurati che il livello di saturazione sia positivo per evitare divisioni per zero o NaN.
    if (saturation_level_ <= 0) {
        saturation_level_ = 1e-6; // Imposta un valore piccolo positivo di default
    }
    // Assicurati che il fattore knee sia positivo per evitare problemi con pow.
    if (knee_factor_ <= 0) {
        knee_factor_ = 1.0; // Imposta un valore di default
    }
}

/**
 * @brief Processa un singolo campione audio attraverso il modello del pentodo.
 *
 * Il modello applica le seguenti fasi di elaborazione:
 * 1.  **Cutoff:** Se il segnale di ingresso è al di sotto della `cutoff_threshold_`,
 * viene limitato a quel valore, simulando il comportamento di cutoff del pentodo.
 * 2.  **Guadagno Lineare:** Il segnale (al di sopra della soglia di cutoff) viene amplificato
 * dal fattore `gain_`.
 * 3.  **Saturazione con "Knee" Effect:** Viene applicata una funzione di saturazione più complessa
 * che combina `tanh` con una potenza, per simulare una compressione più aggressiva e
 * un effetto di "ginocchio" nella curva di distorsione, tipico dei pentodi.
 * 4.  **Offset DC:** Un valore `output_dc_offset_` viene aggiunto all'uscita finale
 * per simulare la tensione di placca a riposo del pentodo.
 *
 * @param input_sample Il campione di ingresso da processare (segnale AC).
 * @return Il campione processato (segnale AC con offset DC).
 */
double Pentode::process(double input_sample) {
    double processed_signal = input_sample;

    // 1. Applica il cutoff
    // Se l'input è sotto la soglia di cutoff, il pentodo non conduce efficacemente.
    processed_signal = std::max(cutoff_threshold_, processed_signal);

    // Calcola il segnale effettivo al di sopra del cutoff
    double signal_above_cutoff = processed_signal - cutoff_threshold_;

    // 2. Applica il guadagno lineare
    double amplified_signal = signal_above_cutoff * gain_;

    // 3. Saturazione con "Knee" Effect
    // Usiamo una combinazione di tanh e una funzione di potenza per una saturazione più dura
    // e un effetto di "ginocchio" (knee) nella curva.
    double saturated_signal;
    if (saturation_level_ > 0) {
        // Normalizza l'input per la funzione di saturazione
        double normalized_signal = amplified_signal / saturation_level_;

        // Applica la funzione di saturazione con effetto "knee"
        // tanh(x) fornisce una saturazione morbida, ma elevando a una potenza dispari
        // si può rendere la curva più ripida e creare un effetto "knee".
        // Per mantenere la simmetria, usiamo std::copysign per conservare il segno.
        saturated_signal = saturation_level_ * std::copysign(std::pow(std::fabs(std::tanh(normalized_signal)), knee_factor_), normalized_signal);
    } else {
        saturated_signal = amplified_signal; // Nessuna saturazione se saturation_level_ è 0
    }

    // 4. Aggiungi l'offset DC all'uscita
    // Questo simula la tensione di placca a riposo del pentodo.
    double final_output = saturated_signal + output_dc_offset_;

    return final_output;
}
