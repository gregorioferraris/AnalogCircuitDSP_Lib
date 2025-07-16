#include "CurrentSource.h"
#include <iostream>

// Costruttore
CurrentSource::CurrentSource(const std::string& name, const std::string& node_plus,
                             const std::string& node_minus, double initial_current)
    : Component(name, node_plus, node_minus), _current_value(initial_current)
{
    // Il costruttore della classe base Component gestisce l'inizializzazione dei nodi.
}

// Implementazione di getStamps
void CurrentSource::getStamps(Eigen::MatrixXd& stamp_A, Eigen::VectorXd& stamp_B,
                           int num_total_equations, double dt,
                           const Eigen::VectorXd& current_solution_guess,
                           const Eigen::VectorXd& prev_solution,
                           double time)
{
    // Una sorgente di corrente contribuisce solo al vettore RHS (stamp_B).
    // stamp_A viene lasciata invariata (tutti zeri per una sorgente di corrente ideale).
    
    // Assicurati che i vettori e le matrici siano della dimensione corretta.
    // In un sistema reale, queste matrici verrebbero passate per riferimento
    // e sarebbero già pre-allocate dal solutore. Qui le inizializziamo
    // come un controllo di sicurezza, ma tipicamente stamp_A e stamp_B
    // sarebbero già globali (o membri di una classe solutore) e verrebbero
    // solo aggiornate qui.
    if (stamp_A.rows() != num_total_equations || stamp_A.cols() != num_total_equations) {
        stamp_A = Eigen::MatrixXd::Zero(num_total_equations, num_total_equations);
    } else {
        stamp_A.setZero(); // Resetta a zero per ogni chiamata (se non lo fa il solutore)
    }

    if (stamp_B.size() != num_total_equations) {
        stamp_B = Eigen::VectorXd::Zero(num_total_equations);
    } else {
        stamp_B.setZero(); // Resetta a zero per ogni chiamata (se non lo fa il solutore)
    }

    // Ottieni gli ID numerici dei nodi
    // Questo è fondamentale per inserire i contributi nelle posizioni corrette della matrice/vettore MNA.
    if (node_ids.size() != 2 || node_ids[0] == -1 || node_ids[1] == -1) {
        // Errore: ID dei nodi non impostati correttamente.
        // Ciò significa che setNodeIds non è stato chiamato o i nodi non sono stati mappati.
        std::cerr << "Errore: ID dei nodi non validi per la sorgente di corrente " << name << std::endl;
        return; 
    }

    int node_plus_id = node_ids[0];
    int node_minus_id = node_ids[1];
    
    double current_val = getCurrent(time); // Ottieni il valore attuale della corrente

    // Applica lo "stamp" della sorgente di corrente al vettore RHS (stamp_B)
    // La corrente ESCE dal nodo positivo (node_plus_id) -> negativo nel RHS
    // La corrente ENTRA nel nodo negativo (node_minus_id) -> positivo nel RHS
    // La riga 0 della matrice MNA corrisponde al nodo di massa (convenzione).
    // Non modifichiamo la riga 0 perché è il riferimento.
    
    if (node_plus_id != 0) { // Se non è il nodo di massa
        stamp_B(node_plus_id) -= current_val;
    }
    if (node_minus_id != 0) { // Se non è il nodo di massa
        stamp_B(node_minus_id) += current_val;
    }
}

// Implementazione di getCurrent (virtuale, può essere sovrascritto)
double CurrentSource::getCurrent(double time) const {
    // Per una sorgente DC o un valore base, restituisce il valore fisso.
    // Classi derivate possono sovrascrivere questo per dipendenza dal tempo (es. sinusoidale, impulso).
    return _current_value;
}

// Implementazione di setCurrent
void CurrentSource::setCurrent(double value) {
    _current_value = value;
}
