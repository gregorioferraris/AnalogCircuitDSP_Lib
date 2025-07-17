// components/Resistor.cpp
#include "Resistor.h"
// #include "../utils/Random.h" // Per applyTolerance - Assumiamo che non sia più necessario o gestito altrove
#include <stdexcept>
#include <iostream> // Per i messaggi di errore

/**
 * @brief Costruttore per il componente Resistor.
 * @param name Il nome univoco del resistore.
 * @param node_names_str Un vettore di stringhe contenente i nomi dei due nodi collegati.
 * @param resistance_nominal Il valore nominale della resistenza in Ohm.
 */
Resistor::Resistor(const std::string& name, const std::vector<std::string>& node_names_str,
                   double resistance_nominal)
    : Component(name, node_names_str) // Chiama il costruttore della classe base
{
    // Non usiamo applyTolerance qui per semplicità, se non è definito o non necessario.
    // resistance = applyTolerance(resistance_nominal, tolerance_percent);
    resistance = resistance_nominal;

    if (resistance <= 0) {
        std::cerr << "Attenzione: Resistenza non positiva per Resistor " << name_ << ". Impostato a 1.0 Ohm." << std::endl;
        resistance = 1.0; // Valore di default positivo
    }
    std::cout << "Resistor " << name_ << " inizializzato con resistenza: " << resistance << " Ohm." << std::endl;
}

/**
 * @brief Applica gli "stamps" del resistore alla matrice MNA (A) e al vettore (B).
 *
 * Per un resistore, gli stamps sono costanti e lineari, basati sulla legge di Ohm.
 *
 * @param num_total_equations Dimensione totale della matrice MNA.
 * @param dt Passo temporale (non usato per i resistori statici).
 * @param x Vettore della soluzione corrente (non usato per i resistori lineari).
 * @param prev_solution La soluzione dal passo temporale precedente (non usata per i resistori statici).
 * @param time Il tempo di simulazione corrente (non usato per i resistori statici).
 * @param A La matrice MNA a cui vengono applicati gli stamps.
 * @param B Il vettore lato destro MNA a cui vengono applicati gli stamps.
 */
void Resistor::getStamps(
    int num_total_equations, double dt,
    const std::vector<double>& x,
    const std::vector<double>& prev_solution,
    double time,
    std::vector<std::vector<double>>& A,
    std::vector<double>& B
) {
    // Ottieni gli indici globali per i due nodi collegati.
    // Assumiamo che node_ids_ contenga [node1, node2]
    if (node_ids_.size() != 2) {
        std::cerr << "Errore: Resistor " << name_ << " si aspetta 2 ID di nodo, ma ne ha " << node_ids_.size() << std::endl;
        return;
    }
    int node1_id = node_ids_[0];
    int node2_id = node_ids_[1];

    // Controlla la validità degli indici
    if (node1_id < 0 || node1_id >= num_total_equations || node2_id < 0 || node2_id >= num_total_equations) {
        std::cerr << "Errore: Indice di nodo non valido per Resistor " << name_ << std::endl;
        return;
    }

    double G = 1.0 / resistance; // Conduttanza

    // Contributi alla matrice MNA (ammettenze)
    // Gli stamps sono:
    // A[node1_id][node1_id] += G
    // A[node2_id][node2_id] += G
    // A[node1_id][node2_id] -= G
    // A[node2_id][node1_id] -= G

    // Applica gli stamps solo se gli indici dei nodi sono validi e diversi da 0 (ground).
    // Il nodo 0 è il ground, e non dovrebbe essere modificato direttamente dagli stamps dei componenti.
    // L'MNA solver dovrebbe già inizializzare A[0][0] = 1 e B[0] = 0.
    // Gli stamps qui sono per i nodi da 1 a N.

    if (node1_id != 0) A[node1_id][node1_id] += G;
    if (node2_id != 0) A[node2_id][node2_id] += G;
    if (node1_id != 0 && node2_id != 0) {
        A[node1_id][node2_id] -= G;
        A[node2_id][node1_id] -= G;
    }
    // Per un resistore, non ci sono contributi al vettore lato destro 'B' dal componente stesso.
}

/**
 * @brief Aggiorna lo stato interno del resistore.
 *
 * Per un resistore ideale, non c'è uno stato interno che evolve nel tempo.
 * Questo metodo è lasciato vuoto per questo modello.
 *
 * @param current_solution Il vettore della soluzione corrente (non usato).
 * @param prev_solution Il vettore della soluzione precedente (non usato).
 * @param dt Il passo temporale (non usato).
 */
void Resistor::updateState(const std::vector<double>& current_solution,
                     const std::vector<double>& prev_solution,
                     double dt) {
    // Nessuno stato interno da aggiornare per un resistore ideale.
    // Questo metodo è intenzionalmente lasciato vuoto.
}
