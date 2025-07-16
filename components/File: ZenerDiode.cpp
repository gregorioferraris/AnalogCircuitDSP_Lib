// File: ZenerDiode.cpp
#include "ZenerDiode.h"
#include <algorithm> // Per std::max (opzionale, per clampare valori)

/**
 * @brief Implementazione del costruttore della classe ZenerDiode.
 * @param name Nome del componente.
 * @param anode_node Nome del nodo dell'anodo.
 * @param cathode_node Nome del nodo del catodo.
 * @param ground_node Nome del nodo di massa.
 * @param zener_voltage La tensione Zener (Vz) in volt.
 * @param forward_voltage La caduta di tensione in polarizzazione diretta (Vf) in volt.
 * @param series_resistance La resistenza in serie (Rs) in ohm quando il diodo conduce.
 * @param reverse_leakage_resistance La resistenza di dispersione inversa (Rr) in ohm quando il diodo è bloccato.
 */
ZenerDiode::ZenerDiode(const std::string& name, const std::string& anode_node, const std::string& cathode_node,
                       const std::string& ground_node, double zener_voltage, double forward_voltage,
                       double series_resistance, double reverse_leakage_resistance)
    : Component(name, anode_node, cathode_node, ground_node), // L'input_node è l'anodo, l'output_node è il catodo
      zener_voltage_(zener_voltage),
      forward_voltage_(forward_voltage),
      series_resistance_(series_resistance),
      reverse_leakage_resistance_(reverse_leakage_resistance)
{
    // Assicurati che i parametri siano validi
    if (zener_voltage_ <= 0) zener_voltage_ = 5.1; // Valore di default ragionevole
    if (forward_voltage_ <= 0) forward_voltage_ = 0.7; // Valore di default per diodo al silicio
    if (series_resistance_ <= 0) series_resistance_ = 1e-3; // Molto piccola ma non zero
    if (reverse_leakage_resistance_ <= 0) reverse_leakage_resistance_ = 1e9; // Molto grande ma non infinita
}

/**
 * @brief Elabora un singolo campione di tensione attraverso il diodo Zener.
 *
 * Questo metodo calcola la corrente che scorre attraverso il diodo Zener
 * in base alla tensione applicata tra anodo e catodo (`voltage_across_diode`).
 *
 * Il modello è piecewise lineare:
 * 1.  **Polarizzazione diretta (V_AK > V_forward):**
 * Il diodo conduce e la corrente è limitata dalla resistenza in serie.
 * Corrente = (V_AK - V_forward) / R_series
 * 2.  **Breakdown Zener (V_AK < -V_zener):**
 * Il diodo conduce in direzione inversa e la tensione ai suoi capi si stabilizza
 * alla tensione Zener. La corrente inversa è limitata dalla resistenza in serie.
 * Corrente = (V_AK + V_zener) / R_series (Nota: V_AK è negativo, quindi la corrente sarà negativa)
 * 3.  **Regione di blocco (tra -V_zener e V_forward):**
 * Il diodo è bloccato. Viene modellato con una resistenza di dispersione molto alta,
 * risultando in una corrente molto piccola.
 * Corrente = V_AK / R_reverse_leakage
 *
 * @param voltage_across_diode La tensione tra anodo e catodo (V_anodo - V_catodo).
 * @return La corrente che scorre attraverso il diodo (A).
 */
double ZenerDiode::process(double voltage_across_diode) {
    double current = 0.0;

    if (voltage_across_diode > forward_voltage_) {
        // Polarizzazione diretta: il diodo conduce
        current = (voltage_across_diode - forward_voltage_) / series_resistance_;
    } else if (voltage_across_diode < -zener_voltage_) {
        // Breakdown Zener: il diodo conduce in direzione inversa
        // La tensione ai capi della resistenza in serie è V_AK - (-Vz) = V_AK + Vz
        current = (voltage_across_diode + zener_voltage_) / series_resistance_;
    } else {
        // Regione di blocco (inversa o piccola diretta): solo corrente di dispersione
        current = voltage_across_diode / reverse_leakage_resistance_;
    }

    return current;
}

// Implementazioni dei metodi setter
void ZenerDiode::setZenerVoltage(double vz) {
    if (vz > 0) zener_voltage_ = vz;
}

void ZenerDiode::setForwardVoltage(double vf) {
    if (vf > 0) forward_voltage_ = vf;
}

void ZenerDiode::setSeriesResistance(double rs) {
    if (rs > 0) series_resistance_ = rs;
}

void ZenerDiode::setReverseLeakageResistance(double rr) {
    if (rr > 0) reverse_leakage_resistance_ = rr;
}
