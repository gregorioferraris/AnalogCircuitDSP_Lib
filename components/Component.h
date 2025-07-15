// components/Component.h
#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>
#include <vector>
#include <numeric>   // Per std::iota (non usato direttamente qui, ma comune per ID)
#include <cmath>     // Per std::pow, std::log10, std::exp, std::tanh
#include <random>    // Per std::mt19937, std::uniform_real_distribution, std::normal_distribution
#include <stdexcept> // Per std::invalid_argument

class Component {
protected:
    std::string name;
    std::vector<std::string> nodeNamesStr; // Nomi dei nodi come stringhe (es. "input", "output")
    std::vector<int> nodeIds;             // ID numerici dei nodi assegnati dal Circuit
    int componentId;                      // ID univoco del componente

    // Generatore di numeri casuali statico per le tolleranze
    static std::mt19937 rng;

public:
    // Costruttore variadico per accettare un numero variabile di nomi di nodi
    template<typename... Args>
    Component(const std::string& name, Args... node_names_str)
        : name(name), componentId(-1) { // -1 indica non ancora assegnato
        (nodeNamesStr.push_back(node_names_str), ...);
    }

    // Distruttore virtuale per polimorfismo
    virtual ~Component() = default;

    // Metodo virtuale puro per clonare il componente (necessario per i sottocircuiti)
    // Ogni classe derivata dovrà implementare questo metodo.
    virtual Component* clone() const = 0;

    // Metodi per impostare gli ID (chiamati da Circuit)
    void setNodeIds(const std::vector<int>& ids) {
        nodeIds = ids;
    }

    void setComponentId(int id) {
        componentId = id;
    }

    // Metodi getter
    const std::string& getName() const { return name; }
    const std::vector<std::string>& getNodeNamesStr() const { return nodeNamesStr; }
    const std::vector<int>& getNodeIds() const { return nodeIds; }
    int getComponentId() const { return componentId; }

    // Metodo per applicare le tolleranze di produzione
    // Utilizza un generatore di numeri casuali statico
    double applyTolerance(double nominal_value, double tolerance_percent, const std::string& distribution = "uniform") {
        if (tolerance_percent == 0.0) {
            return nominal_value;
        }

        if (distribution == "uniform") {
            std::uniform_real_distribution<double> dist(-tolerance_percent / 100.0, tolerance_percent / 100.0);
            double variation_factor = dist(rng);
            return nominal_value * (1.0 + variation_factor);
        } else if (distribution == "normal") {
            // Assumiamo che la tolleranza percentuale rappresenti 3 deviazioni standard (3-sigma)
            double std_dev = nominal_value * (tolerance_percent / (3.0 * 100.0));
            std::normal_distribution<double> dist(nominal_value, std_dev);
            return dist(rng);
        } else {
            throw std::invalid_argument("Tipo di distribuzione non supportato. Usare 'uniform' o 'normal'.");
        }
    }

    // Metodo virtuale puro per ottenere i "contributi" del componente alla matrice MNA e al vettore RHS.
    // Deve essere implementato dalle sottoclassi.
    // current_solution_guess e prev_solution sono vettori di tensioni ai nodi e correnti ausiliarie.
    virtual void getStamps(int num_total_equations, double dt, const std::vector<double>& current_solution_guess, const std::vector<double>& prev_solution, double time,
                           std::vector<std::vector<double>>& A, std::vector<double>& B) = 0;
};

// L'inizializzazione del generatore di numeri casuali statico (Component::rng)
// deve avvenire in un singolo file .cpp, tipicamente in Component.cpp.

#endif // COMPONENT_H
