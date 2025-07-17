// analog_stages/AnalogOutputStage.cpp
#include "AnalogOutputStage.h"
#include "Resistor.h"
#include "Capacitor.h"
#include <iostream>
#include <algorithm> // Per std::find

// Implementazione del metodo addNode
void AnalogOutputStage::addNode(const std::string& node_name) {
    // Aggiunge il nodo solo se non è già presente
    if (std::find(nodes_.begin(), nodes_.end(), node_name) == nodes_.end()) {
        nodes_.push_back(node_name);
        std::cout << "Nodo aggiunto a AnalogOutputStage '" << name_ << "': " << node_name << std::endl;
    }
}

// Costruttore della classe AnalogOutputStage
AnalogOutputStage::AnalogOutputStage(const std::string& name,
                                   const std::string& input_node,
                                   const std::string& output_node,
                                   const std::string& ground_node,
                                   double output_resistance,
                                   double filter_C)
    : name_(name),
      input_nodes_({input_node}),
      output_nodes_({output_node}),
      ground_node_(ground_node)
{
    // Nodi interni specifici per questa istanza del sottocircuito
    filter_node_ = name_ + "_filter_mid"; // Nodo intermedio per il filtro RC

    // Aggiungi i nodi al sottocircuito
    addNode(input_node);
    addNode(output_node);
    addNode(ground_node_);
    addNode(filter_node_);

    std::cout << "Inizializzazione AnalogOutputStage: " << name_ << std::endl;

    // 1. Resistenza di uscita
    // Collega il nodo di ingresso dello stadio di uscita al nodo intermedio del filtro
    // tramite una resistenza che rappresenta l'impedenza di uscita dello stadio precedente
    // o una resistenza di isolamento.
    addComponent<Resistor>(name_ + "_ROutput", {input_node, filter_node_}, output_resistance);
    std::cout << "Aggiunto Resistor " << name_ << "_ROutput" << std::endl;

    // 2. Filtro Passa-Basso (RC) opzionale all'uscita
    // Se filter_C è > 0, aggiungiamo un condensatore in parallelo alla resistenza di uscita
    // per creare un filtro passa-basso. Questo può servire come filtro anti-aliasing
    // o per modellare la risposta in frequenza dell'uscita.
    if (filter_C > 0) {
        addComponent<Capacitor>(name_ + "_COutputFilter", {filter_node_, ground_node_}, filter_C);
        std::cout << "Aggiunto Capacitor " << name_ << "_COutputFilter" << std::endl;
    } else {
        std::cout << "Nessun condensatore di filtro aggiunto (filter_C <= 0)." << std::endl;
    }

    // 3. Collegamento al nodo di uscita esterno
    // Collega il nodo intermedio del filtro al nodo di uscita effettivo del sottocircuito.
    addComponent<Resistor>(name_ + "_RConnectOut", {filter_node_, output_node}, 1e-6); // Quasi un corto
    std::cout << "Aggiunto Resistor " << name_ << "_RConnectOut" << std::endl;

    std::cout << "AnalogOutputStage '" << name_ << "' completamente inizializzato." << std::endl;
}
