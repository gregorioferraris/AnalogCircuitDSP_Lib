// File: ZenerDiode.h
#ifndef ZENER_DIODE_H
#define ZENER_DIODE_H

#include "Component.h" // Presuppone che Component.h definisca la classe base Component
#include <string>

/**
 * @class ZenerDiode
 * @brief Modella un componente diodo Zener semplificato.
 *
 * Questa classe simula il comportamento di un diodo Zener, inclusa la polarizzazione diretta,
 * la regione di blocco inversa e la regione di breakdown Zener.
 * Il metodo process accetta la tensione tra anodo e catodo e restituisce la corrente.
 */
class ZenerDiode : public Component {
public:
    /**
     * @brief Costruttore per la classe ZenerDiode.
     * @param name Nome del componente.
     * @param anode_node Nome del nodo dell'anodo.
     * @param cathode_node Nome del nodo del catodo.
     * @param ground_node Nome del nodo di massa.
     * @param zener_voltage La tensione Zener (Vz) in volt.
     * @param forward_voltage La caduta di tensione in polarizzazione diretta (Vf) in volt.
     * @param series_resistance La resistenza in serie (Rs) in ohm quando il diodo conduce.
     * @param reverse_leakage_resistance La resistenza di dispersione inversa (Rr) in ohm quando il diodo è bloccato.
     */
    ZenerDiode(const std::string& name, const std::string& anode_node, const std::string& cathode_node,
               const std::string& ground_node, double zener_voltage, double forward_voltage = 0.7,
               double series_resistance = 1.0, double reverse_leakage_resistance = 1e9);

    /**
     * @brief Elabora un singolo campione di tensione attraverso il diodo Zener.
     * @param voltage_across_diode La tensione tra anodo e catodo (V_anodo - V_catodo).
     * @return La corrente che scorre attraverso il diodo (positiva dall'anodo al catodo).
     */
    double process(double voltage_across_diode) override;

    // Metodi setter per i parametri
    void setZenerVoltage(double vz);
    void setForwardVoltage(double vf);
    void setSeriesResistance(double rs);
    void setReverseLeakageResistance(double rr);

    // Metodi getter
    double getZenerVoltage() const { return zener_voltage_; }
    double getForwardVoltage() const { return forward_voltage_; }
    double getSeriesResistance() const { return series_resistance_; }
    double getReverseLeakageResistance() const { return reverse_leakage_resistance_; }

private:
    double zener_voltage_;              ///< Tensione di breakdown Zener (Vz).
    double forward_voltage_;            ///< Caduta di tensione in polarizzazione diretta (Vf).
    double series_resistance_;          ///< Resistenza in serie quando conduce.
    double reverse_leakage_resistance_; ///< Resistenza di dispersione in polarizzazione inversa (prima del breakdown).
};

#endif // ZENER_DIODE_H
