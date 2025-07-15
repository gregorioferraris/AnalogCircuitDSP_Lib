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

// Forward declaration per Subcircuit se decidi di usarlo in C++
// class Subcircuit;

class Circuit {
private:
    std::map<std::string, int> nodes; // Mappa: nome del nodo -> ID numerico
    int nextNodeId;                   // Prossimo ID numerico disponibile per un nodo
    std::vector<std::shared_ptr<Component>> components; // Lista di tutti i componenti nel circuito

    std::vector<std::shared_ptr<VoltageSource>> voltageSources; // Lista delle sorgenti di tensione
    std::vector<std::shared_ptr<Splitter>> splitters;           // Lista degli splitter

    // Se si usano sottocircuiti, potrebbero essere qui
    // std::vector<std::shared_ptr<Subcircuit>> subcircuits;

public:
    Circuit();

    // Aggiunge un nodo al circuito se non esiste già e restituisce il suo ID.
    int addNode(const std::string& node_name);

    // Aggiunge un componente al circuito e assegna gli ID numerici ai suoi nodi.
    void addComponent(std::shared_ptr<Component> component);

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
        // Se il nodo non esiste, potresti voler aggiungere un errore o un default
        std::cerr << "Errore: Nodo '" << node_name << "' non trovato nel circuito." << std::endl;
        return -1; // O lancia un'eccezione
    }

    // Metodo per connettere blocchi (se si implementano sottocircuiti in C++)
    // Per ora, questo è un placeholder o una nota concettuale.
    // La logica di connessione blocchi è complessa in C++ e spesso gestita a livello di MnaSolver
    // o tramite un builder di circuito.
    // void connectBlocks(std::shared_ptr<Subcircuit> source_block, std::shared_ptr<Subcircuit> dest_block,
    //                    const std::vector<std::string>& source_output_names, const std::vector<std::string>& dest_input_names);
};

#endif // CIRCUIT_H
