// circuits/DiodeBridgeCompressorSidechain.cpp
#include "DiodeBridgeCompressorSidechain.h"
#include <iostream>
#include <cmath> // Per M_PI

DiodeBridgeCompressorSidechain::DiodeBridgeCompressorSidechain(
    const std::string& name,
    const std::string& input_node,
    const std::string& output_node,
    const std::string& ground_node,
    double diode_is,
    double diode_n,
    double attack_ms,
    double release_ms,
    double filter_cap_ref,
    double sample_rate)
    : Subcircuit(name, input_node, output_node, ground_node), // Le porte esterne sono input_node, output_node, ground_node
      sample_rate_(sample_rate),
      diode_is_(diode_is),
      diode_n_(diode_n),
      filter_C_val_(filter_cap_ref) // Usa la capacità di riferimento
{
    std::cout << "Inizializzazione DiodeBridgeCompressorSidechain: " << name << std::endl;

    // Calcola la resistenza equivalente per il filtro RC.
    // Per un semplice rilevatore di inviluppo RC, attacco e rilascio sono spesso modellati
    // dalla stessa costante di tempo R*C. Se distinti, implica un circuito più complesso
    // (es. percorsi di carica/scarica diversi).
    // Usiamo il tempo di rilascio per la resistenza principale del filtro, poiché è solitamente più lungo.
    // Costante di tempo tau = R*C. R = tau / C. tau in secondi.
    double tau_release_s = release_ms / 1000.0;
    filter_R_val_ = tau_release_s / filter_C_val_;

    // Assicurati che R e C siano positivi e ragionevoli
    if (filter_R_val_ <= 0) filter_R_val_ = 1.0; // Valore minimo di 1 Ohm
    if (filter_C_val_ <= 0) filter_C_val_ = 1e-9; // Valore minimo di 1 nF

    _add_nodes();
    _add_components();
    _connect_nodes();
}

void DiodeBridgeCompressorSidechain::_add_nodes() {
    // Nodi interni per il ponte di diodi
    db_pos_node_ = getName() + "_DB_Pos";
    db_neg_node_ = getName() + "_DB_Neg";
    filter_in_node_ = getName() + "_Filter_In";
    // Il nodo di output esterno è già definito dalla classe base Subcircuit
}

void DiodeBridgeCompressorSidechain::_add_components() {
    // Ponte di Diodi Rettificatore (a onda intera)
    // D1: Ingresso -> DB_Pos
    addInternalComponent(std::make_shared<Diode>(
        getName() + "_D1", getPortNames()[0], db_pos_node_, getPortNames()[2], diode_is_, diode_n_));
    // D2: Ingresso -> DB_Neg
    addInternalComponent(std::make_shared<Diode>(
        getName() + "_D2", getPortNames()[0], db_neg_node_, getPortNames()[2], diode_is_, diode_n_));
    // D3: DB_Pos -> Ingresso
    addInternalComponent(std::make_shared<Diode>(
        getName() + "_D3", db_pos_node_, getPortNames()[0], getPortNames()[2], diode_is_, diode_n_));
    // D4: DB_Neg -> Ingresso
    addInternalComponent(std::make_shared<Diode>(
        getName() + "_D4", db_neg_node_, getPortNames()[0], getPortNames()[2], diode_is_, diode_n_));


    // Filtro RC per il rilevamento dell'inviluppo
    // Resistore dall'uscita rettificata all'ingresso del filtro
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_R_Filter", db_pos_node_, filter_in_node_, filter_R_val_));
    // Condensatore dall'ingresso del filtro a massa
    addInternalComponent(std::make_shared<Capacitor>(
        getName() + "_C_Filter", filter_in_node_, getPortNames()[2], filter_C_val_));

    // L'uscita della side-chain è la tensione attraverso C_filter.
    // Questo è gestito implicitamente mappando il nodo di output esterno al nodo interno filter_in_node_
    // tramite il meccanismo delle porte del Subcircuit.
}

void DiodeBridgeCompressorSidechain::_connect_nodes() {
    // Connetti il nodo di output del filtro RC alla porta di output esterna del sottocircuito.
    // getPortNames()[1] è il nodo di output esterno.
    connectInternalNodes(filter_in_node_, getPortNames()[1]);

    // Connetti il nodo negativo del ponte di diodi a massa (se non già implicitamente gestito)
    connectInternalNodes(db_neg_node_, getPortNames()[2]);
}
