#ifndef CURRENT_SOURCE_H
#define CURRENT_SOURCE_H

#include <string>
#include <vector>
#include <map> // Per mappare i nomi dei nodi agli ID
#include <memory> // Per std::unique_ptr se lo usi nella gestione dei componenti

// Includi Eigen per le operazioni su matrici e vettori
// Assicurati che Eigen sia installato e configurato nel tuo progetto.
// Puoi scaricarlo da https://eigen.tuxfamily.org/
#include <Eigen/Dense>

// --- Classe Base Componente (Versione più funzionale per MNA) ---
// Questa è la classe Component di base da cui CurrentSource erediterà.
// Presuppone un sistema che assegna ID numerici ai nodi per la MNA.
class Component {
public:
    std::string name;
    // Nodi del componente (ad esempio, node_plus, node_minus)
    std::vector<std::string> node_names;
    // ID numerici dei nodi, assegnati da un gestore di nodi esterno
    std::vector<int> node_ids; // Assumiamo che 0 sia il nodo di massa

    Component(const std::string& name, const std::string& node1, const std::string& node2)
        : name(name), node_names({node1, node2}) {}
    
    // Costruttore per componenti con più o meno di 2 nodi (se necessario in futuro)
    Component(const std::string& name, const std::vector<std::string>& nodes)
        : name(name), node_names(nodes) {}

    virtual ~Component() = default;

    // Metodo per impostare gli ID dei nodi una volta che sono stati mappati dal solutore
    void setNodeIds(const std::map<std::string, int>& node_name_to_id_map) {
        node_ids.clear();
        for (const auto& name : node_names) {
            auto it = node_name_to_id_map.find(name);
            if (it != node_name_to_id_map.end()) {
                node_ids.push_back(it->second);
            } else {
                // Gestire l'errore: nodo non trovato nella mappa
                // Per ora, useremo un placeholder invalido o lanceremo un'eccezione
                node_ids.push_back(-1); // ID invalido
                // throw std::runtime_error("Node '" + name + "' not found in map for component " + this->name);
            }
        }
    }

    // Metodi virtuali puri per ottenere i "contributi" del componente
    // al sistema di equazioni MNA.
    // num_total_equations: dimensione totale della matrice e del vettore MNA
    // dt: passo temporale per la simulazione transitoria (per condensatori, induttori)
    // current_solution_guess: soluzione attuale (per componenti non lineari)
    // prev_solution: soluzione dal passo temporale precedente
    // time: tempo corrente della simulazione
    virtual void getStamps(Eigen::MatrixXd& stamp_A, Eigen::VectorXd& stamp_B,
                           int num_total_equations, double dt,
                           const Eigen::VectorXd& current_solution_guess,
                           const Eigen::VectorXd& prev_solution,
                           double time) = 0;
};

// --- Classe CurrentSource ---
class CurrentSource : public Component {
public:
    // Costruttore
    CurrentSource(const std::string& name, const std::string& node_plus,
                  const std::string& node_minus, double initial_current = 0.0);

    // Implementazione del metodo virtuale getStamps dalla classe base Component
    void getStamps(Eigen::MatrixXd& stamp_A, Eigen::VectorXd& stamp_B,
                   int num_total_equations, double dt,
                   const Eigen::VectorXd& current_solution_guess,
                   const Eigen::VectorXd& prev_solution,
                   double time) override;

    // Metodo per ottenere il valore della corrente (può essere sovrascritto)
    virtual double getCurrent(double time) const;

    // Metodo per impostare il valore della corrente
    void setCurrent(double value);

protected:
    double _current_value; // Il valore della corrente (potrebbe essere DC o base AC)
};

#endif // CURRENT_SOURCE_H
