// circuit_solver/Circuit.h
#ifndef CIRCUIT_H
#define CIRCUIT_H

#include <string>
#include <vector>
#include <map>
#include <memory> // Per std::shared_ptr
#include <iostream> // Per debug

#include "components/Component.h"
#include "components/VoltageSource.h" // Per identificare le sorgenti di tensione
#include "components/Splitter.h"      // Per identificare gli splitter

// Forward declaration di Subcircuit (ora è una classe a sé stante)
class Subcircuit; // Dichiarazione forward

class Circuit {
private:
    std::map<std::string, int> nodes; // Mappa: nome del nodo -> ID numerico
    int nextNodeId;                   // Prossimo ID numerico disponibile per un nodo
    std::vector<std::shared_ptr<Component>> components; // Lista di tutti i componenti nel circuito

    std::vector<std::shared_ptr<VoltageSource>> voltageSources; // Lista delle sorgenti di tensione
    std::vector<std::shared_ptr<Splitter>> splitters;           // Lista degli splitter

    // Se si usano sottocircuiti, potrebbero essere qui
    // std::vector<std::shared_ptr<Subcircuit>> subcircuits; // Non direttamente, ma le loro istanze appiattite

public:
    Circuit();

    // Aggiunge un nodo al circuito se non esiste già e restituisce il suo ID.
    int addNode(const std::string& node_name);

    // Aggiunge un componente al circuito e assegna gli ID numerici ai suoi nodi.
    void addComponent(std::shared_ptr<Component> component);

    // NUOVO METODO: Aggiunge un'istanza di un sottocircuito al circuito principale.
    // subcircuit_def: La definizione del sottocircuito da istanziare.
    // instance_connections: Una mappa che collega i nomi delle porte del sottocircuito
    //                       ai nomi dei nodi reali nel circuito principale.
    //                       Esempio: {"input_port": "node_A", "output_port": "node_B"}
    void addSubcircuitInstance(const Subcircuit& subcircuit_def,
                               const std::map<std::string, std::string>& instance_connections);

    // Metodi getter
    int getNumNodes() const { return nextNodeId; } // Il numero di nodi è il prossimo ID disponibile
    const std::map<std::string, int>& getNodeMap() const { return nodes; }
    const std::vector<std::shared_ptr<Component>>& getComponents() const { return components; }
    const std::vector<std::shared_ptr<VoltageSource>>& getVoltageSources() const { return voltageSources; }
    const std::vector<std::shared_ptr<Splitter>>& getSplitters() const { return splitters; }

    // Metodo per ottenere l'ID di un nodo dato il suo nome
    int getNodeId(const std::string& node_name) const {
        auto it = nodes.find(node_name);
        if (it != nodes.end()) {
            return it->second;
        }
        std::cerr << "Errore: Nodo '" << node_name << "' non trovato nel circuito." << std::endl;
        return -1; // O lancia un'eccezione
    }
};

#endif // CIRCUIT_H
