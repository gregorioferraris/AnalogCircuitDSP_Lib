// File: Triode.h
#ifndef TRIODE_H
#define TRIODE_H

#include "Component.h" // Assumiamo che Component.h definisca la classe base Component
#include <string>
#include <cmath>     // Per funzioni matematiche come atan

/**
 * @class Triode
 * @brief Modella un effetto di amplificazione a triodo semplificato con guadagno, saturazione e cutoff.
 *
 * Questa classe simula le caratteristiche non lineari di un triodo,
 * inclusa l'amplificazione, la saturazione del segnale e il cutoff.
 * Non è una simulazione circuitale completa, ma un modello di effetto audio.
 */
class Triode : public Component {
public:
    /**
     * @brief Costruttore per la classe Triode.
     * @param name Nome del componente.
     * @param input_node Nome del nodo di ingresso.
     * @param output_node Nome del nodo di uscita.
     * @param ground_node Nome del nodo di massa.
     * @param gain Fattore di guadagno lineare del triodo.
     * @param saturation_threshold Livello di ingresso a cui inizia la saturazione.
     * @param cutoff_threshold Livello di ingresso sotto il quale il segnale viene tagliato (cutoff).
     * @param output_dc_offset Offset DC aggiunto all'uscita, per simulare la tensione di placca a riposo.
     */
    Triode(const std::string& name, const std::string& input_node, const std::string& output_node,
           const std::string& ground_node, double gain = 10.0, double saturation_threshold = 1.0,
           double cutoff_threshold = -0.5, double output_dc_offset = 50.0);

    /**
     * @brief Processa un singolo campione audio attraverso il modello del triodo.
     * @param input_sample Il campione di ingresso da processare (rappresenta il segnale AC alla griglia).
     * @return Il campione processato (rappresenta il segnale AC alla placca, con offset DC).
     */
    double process(double input_sample) override;

    // Metodi getter per accedere ai parametri (opzionale, ma buona pratica)
    double getGain() const { return gain_; }
    double getSaturationThreshold() const { return saturation_threshold_; }
    double getCutoffThreshold() const { return cutoff_threshold_; }
    double getOutputDcOffset() const { return output_dc_offset_; }

private:
    double gain_;                 ///< Guadagno lineare del triodo.
    double saturation_threshold_; ///< Soglia di saturazione del segnale.
    double cutoff_threshold_;     ///< Soglia di cutoff del segnale di ingresso.
    double output_dc_offset_;     ///< Offset DC aggiunto all'uscita.
};

#endif // TRIODE_H
