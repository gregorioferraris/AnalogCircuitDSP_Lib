// circuits/LowShelfCircuit.cpp
#include "LowShelfCircuit.h"
#include <iostream>
#include <cmath> // Per M_PI e log10

LowShelfCircuit::LowShelfCircuit(const std::string& name,
                                 const std::string& input_node,
                                 const std::string& output_node,
                                 const std::string& ground_node,
                                 double R1_val,
                                 double R2_val,
                                 double C1_val)
    : Subcircuit(name, input_node, output_node, ground_node), // Le porte esterne sono input, output, ground
      R1_val_(R1_val), R2_val_(R2_val), C1_val_(C1_val)
{
    std::cout << "Inizializzazione LowShelfCircuit: " << name << std::endl;

    // Definisci i nomi dei nodi interni
    opamp_vout_node_ = getName() + "_OpAmp_Vout";
    opamp_vminus_node_ = getName() + "_OpAmp_Vminus";
    opamp_vplus_node_ = getName() + "_OpAmp_Vplus";

    // --- Aggiungi i componenti interni al circuito interno del sottocircuito ---

    // Op-Amp
    addInternalComponent(std::make_shared<OpAmp>(
        getName() + "_OpAmp", opamp_vplus_node_, opamp_vminus_node_, opamp_vout_node_));

    // R1 in ingresso
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_R1", getPortNames()[0], opamp_vminus_node_, R1_val_));

    // C1 e R2 in parallelo nel feedback
    // C1 tra OpAmp_Vout e OpAmp_Vminus
    addInternalComponent(std::make_shared<Capacitor>(
        getName() + "_C1", opamp_vout_node_, opamp_vminus_node_, C1_val_));
    // R2 tra OpAmp_Vout e OpAmp_Vminus
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_R2", opamp_vout_node_, opamp_vminus_node_, R2_val_));

    // Connetti l'ingresso non invertente dell'Op-Amp a Massa
    internalCircuit.connectNodes(opamp_vplus_node_, getPortNames()[2]);

    // Mappa il nodo di uscita interno dell'Op-Amp alla porta di uscita esterna del sottocircuito
    internalCircuit.connectNodes(opamp_vout_node_, getPortNames()[1]);
}

/*
// Esempio di come calculate_parameters potrebbe essere implementato se necessario
double LowShelfCircuit::calculateCutoffFrequency() const {
    if (R2_val_ * C1_val_ == 0) return 0.0;
    return 1.0 / (2.0 * M_PI * R2_val_ * C1_val_);
}

double LowShelfCircuit::calculateGainDb() const {
    if (R1_val_ == 0) return std::numeric_limits<double>::infinity();
    return 20 * std::log10(std::abs(-R2_val_ / R1_val_));
}

void LowShelfCircuit::setParameters(double cutoff_freq, double gain_db) {
    // Simile a HighShelf, aggiorna i valori R/C interni e potenzialmente ricalcola
    std::cerr << "Attenzione: setParameters per LowShelfCircuit non completamente implementato per aggiornamenti MNA dinamici." << std::endl;
    // Esempio di sintesi (semplificato):
    // double C_ref = 1.0e-8; // Condensatore di riferimento fisso
    // R2_val_ = 1.0 / (2.0 * M_PI * cutoff_freq * C_ref);
    // R1_val_ = R2_val_ / std::abs(std::pow(10.0, gain_db / 20.0));
    // C1_val_ = C_ref;
}
*/
