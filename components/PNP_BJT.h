// File: PNP_BJT.h
#ifndef PNP_BJT_H
#define PNP_BJT_H

#include "Component.h" // Si assume che Component.h definisca la classe base Component
#include <string>
#include <cmath>     // Per std::exp

/**
 * @class PNP_BJT
 * @brief Modella un effetto semplificato di Transistor a Giunzione Bipolare (BJT) PNP.
 *
 * Questa classe simula le caratteristiche non lineari corrente-tensione di un BJT PNP,
 * concentrandosi principalmente sul suo comportamento nella regione attiva per l'elaborazione audio.
 * L'input è interpretato come la tensione Base-Emettitore (Vbe), e l'output
 * è la corrente di Collettore (Ic) simulata.
 */
class PNP_BJT : public Component {
public:
    /**
     * @brief Costruttore per la classe PNP_BJT.
     * @param name Nome del componente.
     * @param input_node Nome del nodo di input (es. "tensione_base").
     * @param output_node Nome del nodo di output (es. "corrente_collettore").
     * @param ground_node Nome del nodo di massa o emettitore comune.
     * @param beta Guadagno di corrente (hFE) del transistor.
     * @param saturation_current Corrente di saturazione inversa (Is), tipicamente molto piccola.
     * @param thermal_voltage Tensione termica (kT/q), circa 0.026V a temperatura ambiente.
     * @param ideality_factor Fattore di idealità (n), tipicamente tra 1 e 2.
     */
    PNP_BJT(const std::string& name, const std::string& input_node, const std::string& output_node,
            const std::string& ground_node, double beta = 100.0,
            double saturation_current = 1e-14, double thermal_voltage = 0.026, double ideality_factor = 1.0);

    /**
     * @brief Elabora un singolo campione audio attraverso il modello BJT PNP.
     * @param input_sample Il campione di input, interpretato come la tensione Base-Emettitore (Vbe).
     * Per i PNP, Vbe è tipicamente negativa per la conduzione (V_base < V_emitter).
     * @return La corrente di Collettore (Ic) simulata.
     */
    double process(double input_sample) override;

    // Metodi getter
    double getBeta() const { return beta_; }
    double getSaturationCurrent() const { return saturation_current_; }
    double getThermalVoltage() const { return thermal_voltage_; }
    double getIdealityFactor() const { return ideality_factor_; }

private:
    double beta_;              ///< Guadagno di corrente (hFE).
    double saturation_current_; ///< Corrente di saturazione inversa (Is).
    double thermal_voltage_;   ///< Tensione termica (Vt).
    double ideality_factor_;   ///< Fattore di idealità (n).
};

#endif // PNP_BJT_H
