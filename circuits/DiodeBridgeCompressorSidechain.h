// circuits/DiodeBridgeCompressorSidechain.h
#ifndef DIODE_BRIDGE_COMPRESSOR_SIDECHAIN_H
#define DIODE_BRIDGE_COMPRESSOR_SIDECHAIN_H

#include "circuit_solver/subcircuits/Subcircuit.h"
#include "components/Diode.h"
#include "components/Resistor.h"
#include "components/Capacitor.h"

/**
 * @class DiodeBridgeCompressorSidechain
 * @brief Modella una side-chain semplificata per un compressore, basata su un ponte di diodi e un filtro RC.
 *
 * Questo sottocircuito prende un segnale audio in ingresso, lo rettifica tramite un ponte di diodi
 * e lo livella con un filtro RC per generare una tensione di controllo DC.
 * I tempi di attacco e rilascio sono approssimati dalla costante di tempo del filtro RC.
 */
class DiodeBridgeCompressorSidechain : public Subcircuit {
public:
    /**
     * @brief Costruttore per la side-chain del compressore.
     * @param name Nome univoco dell'istanza della side-chain.
     * @param input_node Nome del nodo di ingresso esterno (segnale audio).
     * @param output_node Nome del nodo di uscita esterno (tensione di controllo DC).
     * @param ground_node Nome del nodo di massa esterno.
     * @param diode_is Corrente di saturazione dei diodi (Is).
     * @param diode_n Fattore di idealità dei diodi (n).
     * @param attack_ms Tempo di attacco desiderato in millisecondi.
     * @param release_ms Tempo di rilascio desiderato in millisecondi.
     * @param filter_cap_ref Valore di riferimento del condensatore per il filtro RC.
     * @param sample_rate Frequenza di campionamento audio in Hz (usata per calcoli di tempo).
     */
    DiodeBridgeCompressorSidechain(const std::string& name,
                                   const std::string& input_node,
                                   const std::string& output_node, // This will be the control voltage output
                                   const std::string& ground_node,
                                   double diode_is,
                                   double diode_n,
                                   double attack_ms,
                                   double release_ms,
                                   double filter_cap_ref,
                                   double sample_rate);

private:
    // Nodi interni per il ponte di diodi e il filtro
    std::string db_pos_node_; // Uscita positiva del ponte di diodi
    std::string db_neg_node_; // Uscita negativa del ponte di diodi (solitamente a massa)
    std::string filter_in_node_; // Ingresso al filtro RC

    double sample_rate_;
    double diode_is_;
    double diode_n_;
    double filter_R_val_; // Resistenza calcolata per il filtro
    double filter_C_val_; // Capacità calcolata per il filtro

    /**
     * @brief Aggiunge i nodi interni specifici per questo sottocircuito.
     */
    void _add_nodes();

    /**
     * @brief Aggiunge i componenti interni (diodi, resistori, condensatori).
     */
    void _add_components();

    /**
     * @brief Connette i nodi interni e le porte esterne.
     */
    void _connect_nodes();
};

#endif // DIODE_BRIDGE_COMPRESSOR_SIDECHAIN_H
