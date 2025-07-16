// File: Pentode.h
#ifndef PENTODE_H
#define PENTODE_H

#include "Component.h" // Assumiamo che Component.h definisca la classe base Component
#include <string>
#include <cmath>     // Per funzioni matematiche come tanh, pow

/**
 * @class Pentode
 * @brief Modella un effetto di amplificazione a pentodo semplificato con guadagno, saturazione e cutoff.
 *
 * Questa classe simula le caratteristiche non lineari di un pentodo,
 * inclusa l'amplificazione, una saturazione più pronunciata e il cutoff.
 * Non è una simulazione circuitale completa, ma un modello di effetto audio.
 */
class Pentode : public Component {
public:
    /**
     * @brief Costruttore per la classe Pentode.
     * @param name Nome del componente.
     * @param input_node Nome del nodo di ingresso.
     * @param output_node Nome del nodo di uscita.
     * @param ground_node Nome del nodo di massa.
     * @param gain Fattore di guadagno lineare del pentodo.
     * @param saturation_level Livello di ingresso a cui inizia la saturazione (più pronunciata).
     * @param cutoff_threshold Livello di ingresso sotto il quale il segnale viene tagliato (cutoff).
     * @param output_dc_offset Offset DC aggiunto all'uscita, per simulare la tensione di placca a riposo.
     * @param knee_factor Fattore che influenza la "durezza" della curva di saturazione (effetto ginocchio).
     */
    Pentode(const std::string& name, const std::string& input_node, const std::string& output_node,
            const std::string& ground_node, double gain = 20.0, double saturation_level = 0.8,
            double cutoff_threshold = -0.3, double output_dc_offset = 100.0,
            double knee_factor = 2.0);

    /**
     * @brief Processa un singolo campione audio attraverso il modello del pentodo.
     * @param input_sample Il campione di ingresso da processare (rappresenta il segnale AC alla griglia).
     * @return Il campione processato (rappresenta il segnale AC alla placca, con offset DC).
     */
    double process(double input_sample) override;

    // Metodi getter per accedere ai parametri (opzionale, ma buona pratica)
    double getGain() const { return gain_; }
    double getSaturationLevel() const { return saturation_level_; }
    double getCutoffThreshold() const { return cutoff_threshold_; }
    double getOutputDcOffset() const { return output_dc_offset_; }
    double getKneeFactor() const { return knee_factor_; }

private:
    double gain_;                 ///< Guadagno lineare del pentodo.
    double saturation_level_;     ///< Soglia/livello di saturazione.
    double cutoff_threshold_;     ///< Soglia di cutoff del segnale di ingresso.
    double output_dc_offset_;     ///< Offset DC aggiunto all'uscita.
    double knee_factor_;          ///< Fattore per l'effetto "ginocchio" nella saturazione.
};

#endif // PENTODE_H
