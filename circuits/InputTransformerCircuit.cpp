// circuits/InputTransformerCircuit.cpp
#include "InputTransformerCircuit.h"
#include <iostream>

InputTransformerCircuit::InputTransformerCircuit(
    const std::string& name,
    const std::string& primary_input_node,
    const std::string& primary_ground_node,
    const std::string& secondary_output_node,
    const std::string& secondary_ground_node,
    double turns_ratio,
    double primary_resistance_ohm,
    double secondary_resistance_ohm)
    // Le porte esterne sono input_primario, output_secondario, massa_primario, massa_secondario
    : Subcircuit(name, primary_input_node, secondary_output_node, primary_ground_node, secondary_ground_node),
      turns_ratio_(turns_ratio),
      primary_resistance_ohm_(primary_resistance_ohm),
      secondary_resistance_ohm_(secondary_resistance_ohm)
{
    std::cout << "Inizializzazione InputTransformerCircuit: " << name << std::endl;
    _add_nodes();
    _add_components();
    _connect_nodes();
}

void InputTransformerCircuit::_add_nodes() {
    primary_res_out_node_ = getName() + "_Primary_Res_Out";
    secondary_res_out_node_ = getName() + "_Secondary_Res_Out";
    // I quattro nodi esterni (primary_input, primary_ground, secondary_output, secondary_ground)
    // sono già gestiti come porte dalla base Subcircuit.
}

void InputTransformerCircuit::_add_components() {
    // Aggiungi la resistenza in serie sul primario
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_R_Primary", getPortNames()[0], primary_res_out_node_, primary_resistance_ohm_));

    // Aggiungi la resistenza in serie sul secondario
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_R_Secondary", secondary_res_out_node_, getPortNames()[1], secondary_resistance_ohm_)); // Si connette all'output secondario esterno

    // Aggiungi il componente Trasformatore ideale
    transformer_ = std::make_shared<Transformer>(
        getName() + "_Ideal_Transformer",
        primary_res_out_node_, getPortNames()[2], // Lato primario (dopo R in serie, verso massa primaria)
        secondary_res_out_node_, getPortNames()[3], // Lato secondario (dopo R in serie, verso massa secondaria)
        turns_ratio_);
    addInternalComponent(transformer_);
}

void InputTransformerCircuit::_connect_nodes() {
    // Le connessioni sono implicitamente create dal modo in cui i componenti vengono aggiunti
    // e da come sono definite le porte del Subcircuit.
    // getPortNames()[0] è primary_input_node
    // getPortNames()[1] è secondary_output_node
    // getPortNames()[2] è primary_ground_node
    // getPortNames()[3] è secondary_ground_node

    // Il componente Transformer stesso gestisce le connessioni primarie e secondarie.
    // I resistori in serie sono collegati alle porte esterne.
    // Non sono strettamente necessarie chiamate aggiuntive a connectNodes se la mappatura delle porte
    // e le connessioni dei componenti sono definite correttamente.
    // Tuttavia, se i nodi di massa esterni non sono implicitamente il nodo "0",
    // e sono nodi nominati effettivi, devono essere collegati a "0" se sono massa.
    // Assumendo che getPortNames()[2] e getPortNames()[3] siano *già* mappati alla massa globale '0'
    // o siano intesi come masse flottanti rispetto al sottocircuito.
    // Per un uso generale, è più sicuro collegarli esplicitamente alla massa globale '0' se è l'intenzione.
    internalCircuit.connectNodes(getPortNames()[2], "0"); // Massa primaria alla massa globale
    internalCircuit.connectNodes(getPortNames()[3], "0"); // Massa secondaria alla massa globale
}
