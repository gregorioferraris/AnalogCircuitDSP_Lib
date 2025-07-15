// circuit_solver/subcircuits/Subcircuit.h
#ifndef SUBCIRCUIT_H
#define SUBCIRCUIT_H

#include <string>
#include <vector>
#include <map>
#include <memory> // Per std::shared_ptr

// Forward declaration di Circuit per evitare dipendenze circolari
// Circuit ha bisogno di Subcircuit, e Subcircuit ha bisogno di Circuit (per internalCircuit)
// ma Circuit è già definito nel suo .h.
// In questo caso, Subcircuit ha bisogno della definizione completa di Circuit, quindi lo includiamo.
#include "circuit_solver/Circuit.h" // Includi il Circuit completo

// Questa classe rappresenta la DEFINIZIONE di un sottocircuito riutilizzabile.
// Non è un componente in sé, ma un template per creare istanze nel circuito principale.
class Subcircuit {
private:
    std::string name;
    Circuit internalCircuit; // Il circuito interno che definisce il sottocircuito
    std::vector<std::string> portNames; // Nomi dei nodi che fungono da "porte" esterne

public:
    // Costruttore: prende un nome e i nomi delle porte esterne
    template<typename... Args>
    Subcircuit(const std::string& subcircuit_name, Args... external_port_names)
        : name(subcircuit_name) {
        (portNames.push_back(external_port_names), ...);
    }

    // Metodo per aggiungere componenti al circuito interno del sottocircuito
    void addInternalComponent(std::shared_ptr<Component> component) {
        internalCircuit.addComponent(component);
    }

    // Getter per il nome del sottocircuito
    const std::string& getName() const { return name; }

    // Getter per il circuito interno
    const Circuit& getInternalCircuit() const { return internalCircuit; }

    // Getter per i nomi delle porte
    const std::vector<std::string>& getPortNames() const { return portNames; }
};

#endif // SUBCIRCUIT_H
