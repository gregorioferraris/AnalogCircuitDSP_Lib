// circuits/HighShelfCircuit.h
#ifndef HIGHSHELF_CIRCUIT_H
#define HIGHSHELF_CIRCUIT_H

#include "circuit_solver/subcircuits/Subcircuit.h"
#include "components/Resistor.h"
#include "components/Capacitor.h"
#include "components/OpAmp.h"

/**
 * @class HighShelfCircuit
 * @brief Rappresenta un circuito di filtro High-Shelf attivo.
 *
 * Questo sottocircuito implementa un filtro shelf che aumenta o diminuisce
 * il guadagno delle frequenze al di sopra di una frequenza di taglio specificata.
 * È basato su una configurazione di amplificatore operazionale.
 */
class HighShelfCircuit : public Subcircuit {
public:
    /**
     * @brief Costruttore per un filtro High-Shelf attivo.
     * @param name Nome univoco dell'istanza del filtro.
     * @param input_node Nome del nodo di ingresso esterno.
     * @param output_node Nome del nodo di uscita esterno.
     * @param ground_node Nome del nodo di massa esterno.
     * @param R1_val Valore del resistore R1.
     * @param R2_val Valore del resistore R2.
     * @param C1_val Valore del condensatore C1.
     */
    HighShelfCircuit(const std::string& name,
                     const std::string& input_node,
                     const std::string& output_node,
                     const std::string& ground_node,
                     double R1_val,
                     double R2_val,
                     double C1_val);

    // Metodi per calcolare/impostare i parametri (se necessari per sintesi/analisi)
    // Questi interagirebbero tipicamente con i valori dei componenti interni.
    // Per MNA, questi aggiornerebbero semplicemente i valori interni di R/C, e il circuito
    // dovrebbe essere ricalcolato dal solutore.
    // double calculateCutoffFrequency() const;
    // double calculateGainDb() const;
    // void setParameters(double cutoff_freq, double gain_db);

private:
    // Nomi dei nodi interni
    std::string opamp_vout_node_;
    std::string opamp_vminus_node_;
    std::string opamp_vplus_node_;
    std::string feedback_node_; // Tra R1, C1, R2

    // Parametri (memorizzati se necessari per calculate_parameters/set_parameters)
    double R1_val_;
    double R2_val_;
    double C1_val_;
};

#endif // HIGHSHELF_CIRCUIT_H
