// File: Triode.cpp
#include "Triode.h"
#include <algorithm> // Per std::max, se necessario

/**
 * @brief Implementazione del costruttore della classe Triode.
 * @param name Nome del componente.
 * @param input_node Nome del nodo di ingresso.
 * @param output_node Nome del nodo di uscita.
 * @param ground_node Nome del nodo di massa.
 * @param gain Fattore di guadagno lineare.
 * @param saturation_threshold Livello di ingresso a cui inizia la saturazione.
 * @param cutoff_threshold Livello di ingresso sotto il quale il segnale viene tagliato.
 * @param output_dc_offset Offset DC aggiunto all'uscita.
 */
Triode::Triode(const std::string& name, const std::string& input_node, const std::string& output_node,
               const std::string& ground_node, double gain, double saturation_threshold,
               double cutoff_threshold, double output_dc_offset)
    : Component(name, input_node, output_node, ground_node), // Chiama il costruttore della classe base
      gain_(gain),
      saturation_threshold_(saturation_threshold),
      cutoff_threshold_(cutoff_threshold),
      output_dc_offset_(output_dc_offset)
{
    // Assicurati che la soglia di saturazione sia positiva per evitare divisioni per zero o NaN.
    if (saturation_threshold_ <= 0) {
        saturation_threshold_ = 1e-6; // Imposta un valore piccolo positivo di default
    }
}

/**
 * @brief Processa un singolo campione audio attraverso il modello del triodo.
 *
 * Il modello applica le seguenti fasi di elaborazione:
 * 1.  **Cutoff:** Se il segnale di ingresso è al di sotto della `cutoff_threshold_`,
 * viene limitato a quel valore, simulando il comportamento di cutoff del triodo.
 * 2.  **Guadagno Lineare:** Il segnale (al di sopra della soglia di cutoff) viene amplificato
 * dal fattore `gain_`.
 * 3.  **Saturazione:** Viene applicata una funzione `atan` per simulare una saturazione morbida.
 * Questo comprime i picchi del segnale, tipico delle valvole.
 * 4.  **Offset DC:** Un valore `output_dc_offset_` viene aggiunto all'uscita finale
 * per simulare la tensione di placca a riposo del triodo.
 *
 * @param input_sample Il campione di ingresso da processare (segnale AC).
 * @return Il campione processato (segnale AC con offset DC).
 */
double Triode::process(double input_sample) {
    double processed_signal = input_sample;

    // 1. Applica il cutoff
    // Se l'input è sotto la soglia di cutoff, il triodo non conduce efficacemente.
    // Limitiamo l'input a questa soglia.
    processed_signal = std::max(cutoff_threshold_, processed_signal);

    // Calcola il segnale effettivo al di sopra del cutoff
    double signal_above_cutoff = processed_signal - cutoff_threshold_;

    // 2. Applica il guadagno lineare
    double amplified_signal = signal_above_cutoff * gain_;

    // 3. Saturazione (usando atan per una curva di saturazione morbida)
    // Normalizza il segnale rispetto alla soglia di saturazione prima di applicare atan.
    double saturated_signal;
    if (saturation_threshold_ > 0) {
        saturated_signal = saturation_threshold_ * std::atan(amplified_signal / saturation_threshold_);
    } else {
        // Se la soglia è non positiva, non applicare saturazione.
        saturated_signal = amplified_signal;
    }

    // 4. Aggiungi l'offset DC all'uscita
    // Questo simula la tensione di placca a riposo del triodo.
    double final_output = saturated_signal + output_dc_offset_;

    return final_output;
}
