#ifndef ANALOG_INPUT_STAGE_H
#define ANALOG_INPUT_STAGE_H

#include <string>
#include <vector>
#include <memory> // Per std::unique_ptr

// --- Dichiarazioni Forward per i Componenti ---
// Assumiamo che queste classi esistano o saranno definite altrove nella tua libreria.
// Servono per poter usare puntatori a questi tipi senza includere i loro file .h qui.
class Component; // Classe base per tutti i componenti
class Resistor;
class Capacitor;
class Diode;

// --- Classe Base Componente (Esempio semplificato) ---
// Questa è una versione molto basilare. La tua vera classe Component
// avrà più metodi per la simulazione (e.g., getConductance, getCurrent, ecc.)
class Component {
public:
    std::string name;
    // Nodi collegati (potrebbe essere una lista di stringhe o di indici di nodo)
    // Per semplicità, qui useremo stringhe per i nomi dei nodi, ma per una vera simulazione MNA,
    // avresti bisogno di un modo per mappare questi nomi a indici numerici.
    std::vector<std::string> nodes;

    Component(const std::string& name, const std::vector<std::string>& nodes)
        : name(name), nodes(nodes) {}

    virtual ~Component() = default;

    // Metodi virtuali puri che le classi derivate dovranno implementare
    // Questi sono placeholder per come interagirai con un solutore di circuiti.
    // virtual double getConductance() const = 0;
    // virtual double getCurrent() const = 0;
};

// --- Classi dei Componenti Specifici (Esempi semplificati) ---
// Queste sono solo dichiarazioni. Le implementazioni complete dipenderanno
// da come intendi integrare con il tuo solutore di circuiti C++.

class Resistor : public Component {
public:
    double resistance;
    Resistor(const std::string& name, const std::string& node1, const std::string& node2, double r)
        : Component(name, {node1, node2}), resistance(r) {}
};

class Capacitor : public Component {
public:
    double capacitance;
    Capacitor(const std::string& name, const std::string& node1, const std::string& node2, double c)
        : Component(name, {node1, node2}), capacitance(c) {}
};

class Diode : public Component {
public:
    double Is; // Corrente di saturazione
    double N;  // Coefficiente di emissione
    Diode(const std::string& name, const std::string& anode, const std::string& cathode, double is, double n)
        : Component(name, {anode, cathode}), Is(is), N(n) {}
};

// --- Classe AnalogInputStage ---
class AnalogInputStage {
public:
    // Costruttore con parametri di default
    AnalogInputStage(const std::string& name,
                     const std::string& input_node,
                     const std::string& output_node,
                     const std::string& ground_node = "0",
                     double input_resistance = 10e3,
                     double gain_factor = 1.0,
                     double filter_R = 1e3,
                     double filter_C = 100e-9,
                     double diode_saturation_current = 1e-9,
                     double diode_emission_coefficient = 1.0);

    // Metodo per ottenere il nome del sottocircuito
    std::string getName() const { return name_; }

    // Metodo per accedere ai nodi del sottocircuito (se necessario per il solutore)
    const std::vector<std::string>& getNodes() const { return nodes_; }

    // Metodo per accedere ai componenti (per passarli al solutore)
    const std::vector<std::unique_ptr<Component>>& getComponents() const { return components_; }

    // Metodo per ottenere i nodi di ingresso/uscita
    std::string getInputNode() const { return input_nodes_[0]; } // Assumiamo un solo input
    std::string getOutputNode() const { return output_nodes_[0]; } // Assumiamo un solo output

private:
    std::string name_;
    std::vector<std::string> input_nodes_;
    std::vector<std::string> output_nodes_;
    std::vector<std::string> nodes_; // Tutti i nodi, inclusi quelli interni
    std::vector<std::unique_ptr<Component>> components_; // Usiamo unique_ptr per gestire la memoria

    // Nodi interni
    std::string gain_node_;
    std::string dist_node_;
    std::string filter_node_;
    std::string ground_node_; // Il nodo di riferimento a massa
    std::string diode_anode_node_; // Nodo interno per il diodo
    std::string filter_mid_node_; // Nodo interno per il filtro RC
    
    // Metodo helper per aggiungere un nodo (e assicurarsi che sia unico)
    void addNode(const std::string& node_name);
    
    // Metodo helper per aggiungere un componente
    template<typename T, typename... Args>
    void addComponent(Args&&... args) {
        components_.push_back(std::make_unique<T>(std::forward<Args>(args)...));
    }
};

#endif // ANALOG_INPUT_STAGE_H
