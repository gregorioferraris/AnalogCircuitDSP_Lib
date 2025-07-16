#include "AnalogInputStage.h"
#include <iostream>
#include <algorithm> // Per std::find

// Implementazione del metodo addNode
void AnalogInputStage::addNode(const std::string& node_name) {
    // Aggiunge il nodo solo se non è già presente
    if (std::find(nodes_.begin(), nodes_.end(), node_name) == nodes_.end()) {
        nodes_.push_back(node_name);
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

    // 1. Resistenza di Ingresso
    // La resistenza di ingresso è tra input_node e il nodo di output dello stadio di guadagno
    addComponent<Resistor>(name_ + "_Rin", input_node, gain_node_, input_resistance);

    // 2. Stadio di Guadagno (semplice partitore resistivo per attenuazione/guadagno)
    if (gain_factor < 1.0) {
        // Calcola la resistenza in serie per ottenere l'attenuazione desiderata
        // Nota: questa è una semplificazione. In un partitore vero, avresti un'altra resistenza a massa
        // per creare il rapporto. Qui assumiamo che R_in e R_gain formino il partitore con il carico successivo.
        // Questa logica necessiterebbe di un aggiustamento per un partitore standard.
        // Ad esempio: Vout = Vin * (R_parallel / (R_series + R_parallel))
        // La tua logica python era R_gain_series = input_resistance / gain_factor - input_resistance.
        // Questo implica che R_in è la resistenza superiore e R_gain_series è la resistenza inferiore,
        // ma R_gain_series dovrebbe andare a massa per un partitore tradizionale.
        // Per mantenere la logica Python, lo interpretiamo come una resistenza in serie che modifica l'impedenza vista.
        
        // Se vogliamo un'attenuazione, l'idea è che la R_in faccia parte del partitore.
        // Per ottenere l'attenuazione 'gain_factor', consideriamo che l'impedenza di ingresso è data.
        // Se gain_factor < 1, significa che il segnale viene ridotto.
        // Un semplice modo per fare attenuazione passiva con resistori è un partitore.
        // R_upper = input_resistance (già definita)
        // R_lower = ?
        // Vout = Vin * (R_lower / (R_upper + R_lower)) = Vin * gain_factor
        // gain_factor = R_lower / (R_upper + R_lower)
        // R_lower = gain_factor * (R_upper + R_lower)
        // R_lower = gain_factor * R_upper + gain_factor * R_lower
        // R_lower * (1 - gain_factor) = gain_factor * R_upper
        // R_lower = (gain_factor * R_upper) / (1 - gain_factor)
        
        // Questo calcolo è per un partitore resistivo standard dove R_upper è tra Vin e Vout,
        // e R_lower è tra Vout e massa.
        // La tua logica Python era: r_gain_series = input_resistance / gain_factor - input_resistance.
        // Questo sembra voler modellare una resistenza in serie al segnale dopo R_in per arrivare a un certo punto.
        // Manteniamo la logica Python per coerenza.
        double r_gain_series = input_resistance / gain_factor - input_resistance;
        if (r_gain_series < 1e-6) r_gain_series = 1e-6; // Evita valori troppo piccoli o negativi
        addComponent<Resistor>(name_ + "_RGain", gain_node_, dist_node_, r_gain_series);
    } else {
        // Se guadagno >= 1, consideralo come passante diretto (resistenza molto piccola)
        addComponent<Resistor>(name_ + "_RPass", gain_node_, dist_node_, 1e-6); // Quasi un corto
    }

    // 3. Distorsione (diodo soft-clipping)
    // R_dist_series limita la corrente attraverso il diodo.
    addComponent<Resistor>(name_ + "_R_dist_series", dist_node_, diode_anode_node_, 1e3);
    // Il diodo è tra il suo anodo_node e il nodo di massa
    addComponent<Diode>(name_ + "_D1", diode_anode_node_, ground_node_, diode_saturation_current, diode_emission_coefficient);
    // R_dist_parallel è tra il nodo del diodo e il nodo di ingresso del filtro
    addComponent<Resistor>(name_ + "_R_dist_parallel", diode_anode_node_, filter_node_, 10e3);

    // 4. Filtro Anti-Aliasing Passa-Basso (RC)
    // La resistenza del filtro è tra l'uscita dello stadio di distorsione e il nodo intermedio del filtro
    addComponent<Resistor>(name_ + "_FilterR", filter_node_, filter_mid_node_, filter_R);
    // Il condensatore del filtro è tra il nodo intermedio e la massa
    addComponent<Capacitor>(name_ + "_FilterC", filter_mid_node_, ground_node_, filter_C);
    
    // Output del sottocircuito: collega il nodo di uscita del filtro all'output_node dichiarato
    addComponent<Resistor>(name_ + "_ROut", filter_mid_node_, output_node, 1e-6); // R quasi 0 per collegare a output_node

    std::cout << "Stadio di Ingresso Analogico '" << name_ << "' creato. Impedenze, guadagno e distorsione modellati." << std::endl;
}
