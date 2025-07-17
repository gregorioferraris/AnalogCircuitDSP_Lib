// analog_stages/AnalogInputStage.cpp
#include "AnalogInputStage.h"
#include "Resistor.h"
#include "Capacitor.h"
#include "Diode.h"
#include <iostream>
#include <algorithm> // Per std::find

// Implementazione del metodo addNode
void AnalogInputStage::addNode(const std::string& node_name) {
    // Aggiunge il nodo solo se non è già presente
    if (std::find(nodes_.begin(), nodes_.end(), node_name) == nodes_.end()) {
        nodes_.push_back(node_name);
        std::cout << "Nodo aggiunto a AnalogInputStage '" << name_ << "': " << node_name << std::endl;
    }
}

// Costruttore della classe AnalogInputStage
AnalogInputStage::AnalogInputStage(const std::string& name,
                                 const std::string& input_node,
                                 const std::string& output_node,
                                 const std::string& ground_node,
                                 double input_resistance,
                                 double gain_factor,
                                 double filter_R,
                                 double filter_C,
                                 double diode_saturation_current,
                                 double diode_emission_coefficient)
    : name_(name),
      input_nodes_({input_node}),
      output_nodes_({output_node}),
      ground_node_(ground_node)
{
    // Nodi interni specifici per questa istanza del sottocircuito
    gain_node_ = name_ + "_gain_out";
    dist_node_ = name_ + "_dist_out";
    filter_node_ = name_ + "_filter_out";
    diode_anode_node_ = name_ + "_d1_in";
    filter_mid_node_ = name_ + "_filter_mid";

    // Aggiungi i nodi al sottocircuito
    addNode(input_node);
    addNode(output_node);
    addNode(ground_node_);
    addNode(gain_node_);
    addNode(dist_node_);
    addNode(filter_node_);
    addNode(diode_anode_node_);
    addNode(filter_mid_node_);

    std::cout << "Inizializzazione AnalogInputStage: " << name_ << std::endl;

    // 1. Resistenza di ingresso
    // Collega il nodo di ingresso alla massa tramite la resistenza di ingresso.
    // Questo simula l'impedenza di ingresso dello stadio.
    addComponent<Resistor>(name_ + "_Rin", {input_node, ground_node_}, input_resistance);
    std::cout << "Aggiunto Resistor " << name_ << "_Rin" << std::endl;

    // 2. Stadio di Guadagno (semplificato come un resistore in serie)
    // Se il guadagno è significativo (es. > 1), potremmo modellare una resistenza in serie
    // o un VCVS. Per semplicità, se gain_factor è 1, lo consideriamo un corto.
    // Altrimenti, un resistore in serie può rappresentare una caduta di tensione
    // o un divisore di tensione se combinato con altri componenti.
    // Per un vero "guadagno", avremmo bisogno di un componente attivo (es. OpAmp, Transistor).
    // Qui, per simulare un "guadagno" in modo passivo, useremo un resistore
    // che porta il segnale dal nodo di ingresso al nodo di uscita del guadagno.
    // Questo è un placeholder e andrebbe sostituito con un modello attivo per un vero guadagno.
    if (gain_factor > 0 && gain_factor != 1.0) {
        // Questo è un modello molto semplificato. Un vero guadagno richiederebbe un VCVS o un transistor.
        // Qui, simuliamo una "trasformazione" del segnale.
        // Per ora, lo lascio come un resistore passante.
        // TODO: Implementare un vero stadio di guadagno attivo (es. VCVS, Transistor).
        addComponent<Resistor>(name_ + "_RGain", {input_node, gain_node_}, input_resistance / gain_factor); // Esempio
        std::cout << "Aggiunto Resistor " << name_ << "_RGain" << std::endl;
    } else { // Se guadagno è 1 o non positivo, consideralo come passante diretto (resistenza molto piccola)
        addComponent<Resistor>(name_ + "_RPass", {input_node, gain_node_}, 1e-6); // Quasi un corto
        std::cout << "Aggiunto Resistor " << name_ << "_RPass (passante)" << std::endl;
    }

    // 3. Distorsione (diodo soft-clipping)
    // R_dist_series limita la corrente attraverso il diodo.
    addComponent<Resistor>(name_ + "_R_dist_series", {gain_node_, diode_anode_node_}, 1e3);
    std::cout << "Aggiunto Resistor " << name_ << "_R_dist_series" << std::endl;
    // Il diodo è tra il suo anodo_node e il nodo di massa
    addComponent<Diode>(name_ + "_D1", {diode_anode_node_, ground_node_}, diode_saturation_current, diode_emission_coefficient);
    std::cout << "Aggiunto Diode " << name_ << "_D1" << std::endl;
    // R_dist_parallel è tra il nodo del diodo e il nodo di ingresso del filtro
    addComponent<Resistor>(name_ + "_R_dist_parallel", {diode_anode_node_, dist_node_}, 10e3);
    std::cout << "Aggiunto Resistor " << name_ << "_R_dist_parallel" << std::endl;

    // 4. Filtro Anti-Aliasing Passa-Basso (RC)
    // La resistenza del filtro è tra l'uscita dello stadio di distorsione e il nodo intermedio del filtro
    addComponent<Resistor>(name_ + "_FilterR", {dist_node_, filter_mid_node_}, filter_R);
    std::cout << "Aggiunto Resistor " << name_ << "_FilterR" << std::endl;
    // Il condensatore del filtro è tra il nodo intermedio e la massa
    addComponent<Capacitor>(name_ + "_FilterC", {filter_mid_node_, ground_node_}, filter_C);
    std::cout << "Aggiunto Capacitor " << name_ << "_FilterC" << std::endl;
    
    // Output del sottocircuito: collega il nodo di uscita del filtro all'output_node dichiarato
    addComponent<Resistor>(name_ + "_ROut", {filter_mid_node_, output_node}, 1e-6); // Quasi un corto
    std::cout << "Aggiunto Resistor " << name_ << "_ROut" << std::endl;

    std::cout << "AnalogInputStage '" << name_ << "' completamente inizializzato." << std::endl;
}
