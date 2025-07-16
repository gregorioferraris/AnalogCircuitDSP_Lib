// components/Component.h
#ifndef COMPONENT_H
#define COMPONENT_H

#include <string>
#include <vector>
#include <map>
#include <memory> // For std::shared_ptr
#include <Eigen/Dense> // Include Eigen for matrix/vector types in method signatures

class Component {
public:
    // Costruttore: accetta il nome del componente e una lista di nomi di nodi come stringhe
    Component(const std::string& name, const std::vector<std::string>& nodeNames);
    
    // Distruttore virtuale per consentire la corretta deallocazione delle classi derivate
    virtual ~Component() = default;

    // Metodo virtuale puro per clonare il componente (necessario per i sottocircuiti)
    // Ogni classe derivata dovrà implementare questo metodo.
    virtual Component* clone() const = 0;

    // Metodi per impostare gli ID numerici dei nodi (chiamato da Circuit)
    // Per componenti con nodi ordinati (es. Resistor, Capacitor)
    void setNodeIds(const std::vector<int>& ids);
    // Per componenti con nodi nominati (es. Triode, MOSFET)
    void setNodeIds(const std::map<std::string, int>& ids);

    // Metodo per impostare l'ID univoco del componente (chiamato da Circuit)
    void setComponentId(int id);

    // Metodo virtuale puro per ottenere i "contributi" del componente alla matrice MNA e al vettore RHS.
    // Deve essere implementato dalle sottoclassi.
    // 'x_current_guess' e 'prev_solution' sono le soluzioni del solutore.
    // 'time' è il tempo attuale della simulazione. 'dt' è il passo temporale.
    virtual void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) = 0;

    // Metodo virtuale per calcolare la corrente non lineare (usato da MnaSolver::_system_equations)
    // Questo metodo sarà implementato solo dai componenti non lineari.
    // I parametri dipendono dal tipo di componente (es. Vd per diodo, Vgs/Vds per MOSFET).
    virtual double calculateNonlinearCurrent(double param1, double param2) const {
        // Implementazione di default per componenti lineari o non applicabili
        return 0.0;
    }
    
    // Metodo virtuale per l'aggiornamento dello stato (per C, L, ecc.)
    virtual void updateState(double v_curr, double i_curr) {}

    // Metodo virtuale per ottenere il valore di una sorgente (per VoltageSource, CurrentSource)
    virtual double getValue(double time) const { return 0.0; }
    // Metodo virtuale per impostare il valore di una sorgente (per DelayLine che controlla Vs, o per input audio)
    virtual void setValue(double value) {}

    // Metodo virtuale per l'aggiornamento dei componenti funzionali (es. DelayLine)
    virtual double update(double input_sample) { return input_sample; }
    virtual void reset() {} // Per resettare lo stato di componenti funzionali

    // Getters
    const std::string& getName() const { return name; }
    const std::vector<std::string>& getNodeNamesStr() const { return nodeNamesStr; }
    int getComponentId() const { return componentId; }

    // Metodi per identificare il tipo di nodo ID
    bool hasNamedNodes() const { return !nodeIdsMap.empty(); }
    const std::vector<int>& getNodeIdsVector() const { return nodeIdsVector; }
    const std::map<std::string, int>& getNodeIdsMap() const { return nodeIdsMap; }

protected:
    std::string name;
    std::vector<std::string> nodeNamesStr; // Nomi dei nodi come stringhe (es. {"input", "output"})
    int componentId;

    // Gli ID dei nodi possono essere memorizzati come vettore (per nodi ordinati)
    // o come mappa (per nodi nominati come 'drain', 'gate', 'source').
    // Useremo uno dei due, non entrambi contemporaneamente, a seconda del tipo di componente.
    std::vector<int> nodeIdsVector;
    std::map<std::string, int> nodeIdsMap;

    // Indice per la corrente delle sorgenti di tensione (gestito da MnaSolver)
    int current_index = -1;
    friend class MnaSolver; // MnaSolver può accedere a current_index
};

#endif // COMPONENT_H
