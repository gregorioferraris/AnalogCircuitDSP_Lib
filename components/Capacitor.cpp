// components/Capacitor.cpp
#include "Capacitor.h"
// #include "../utils/Random.h" // Per applyTolerance - Assumendo che non sia più necessario o sia gestito altrove
#include <stdexcept>
#include <iostream> // Per messaggi di debug

/**
 * @brief Costruttore per il componente Capacitor.
 * @param name Il nome univoco del condensatore.
 * @param node_names_str Un vettore di stringhe contenente i nomi dei due nodi collegati.
 * @param capacitance_nominal Il valore nominale della capacità in Farad.
 */
Capacitor::Capacitor(const std::string& name, const std::vector<std::string>& node_names_str,
                     double capacitance_nominal)
    : Component(name, node_names_str), v_prev(0.0), i_prev(0.0) // Chiama il costruttore della classe base
{
    // Non stiamo usando applyTolerance qui per semplicità, se non è definito o non necessario.
    // capacitance = applyTolerance(capacitance_nominal, tolerance_percent);
    capacitance = capacitance_nominal;

    if (capacitance <= 0) {
        std::cerr << "Warning: Capacità non positiva per il Condensatore " << name_ << ". Impostato a 1e-12 F." << std::endl;
        capacitance = 1e-12; // Valore predefinito positivo ragionevole (1 pF)
    }
    std::cout << "Condensatore " << name_ << " inizializzato con capacità: " << capacitance << " F." << std::endl;
}

/**
 * @brief Applica gli "stamps" del condensatore alla matrice MNA (A) e al vettore (B).
 *
 * Utilizza il metodo di integrazione trapezoidale per modellare il condensatore
 * nell'analisi transitoria.
 *
 * @param num_total_equations Dimensione totale della matrice MNA.
 * @param dt Passo temporale.
 * @param x Vettore della soluzione corrente (tensioni dei nodi e correnti di ramo).
 * @param prev_solution La soluzione dal passo temporale precedente.
 * @param time Tempo di simulazione corrente.
 * @param A La matrice MNA a cui vengono applicati gli stamps.
 * @param B Il vettore lato destro MNA a cui vengono applicati gli stamps.
 */
void Capacitor::getStamps(
    int num_total_equations, double dt,
    const std::vector<double>& x,
    const std::vector<double>& prev_solution,
    double time,
    std::vector<std::vector<double>>& A,
    std::vector<double>& B
) {
    if (node_ids_.size() != 2) {
        std::cerr << "Errore: Condensatore " << name_ << " si aspetta 2 ID di nodo." << std::endl;
        return;
    }
    int node1_id = node_ids_[0];
    int node2_id = node_ids_[1];

    // Controlla per indici validi
    if (node1_id < 0 || node1_id >= num_total_equations || node2_id < 0 || node2_id >= num_total_equations) {
        std::cerr << "Errore: Indice di nodo non valido per il Condensatore " << name_ << std::endl;
        return;
    }

    // Calcola la conduttanza equivalente (G_eq) per il modello di integrazione trapezoidale
    // I_C(t) = G_eq * V_C(t) + I_srcC
    // dove G_eq = 2C/dt
    // e I_srcC = (2C/dt) * V_C(t-dt) + I_C(t-dt)
    double G_eq = 2.0 * capacitance / dt;

    // Recupera la tensione e la corrente del condensatore dal passo precedente
    // prev_solution contiene le tensioni dei nodi e le correnti di ramo.
    // Dobbiamo estrarre la tensione tra node1 e node2 e la corrente del condensatore.
    // Assumiamo che v_prev e i_prev siano stati aggiornati in updateState al passo precedente.
    // V_C_prev = prev_solution[node1_id] - prev_solution[node2_id]; // Tensione ai capi del condensatore al passo precedente
    // La corrente del condensatore al passo precedente non è direttamente nel vettore prev_solution
    // a meno che non sia una variabile ausiliaria, ma è memorizzata in i_prev.

    double V_C_prev = v_prev; // Tensione ai capi del condensatore al passo precedente
    double I_C_prev = i_prev; // Corrente attraverso il condensatore al passo precedente

    // Calcola la sorgente di corrente equivalente per il lato destro (B)
    double I_srcC = G_eq * V_C_prev + I_C_prev;

    // Contributi alla matrice MNA (parte dipendente dalla tensione corrente)
    // Gli stamps sono:
    // A[node1_id][node1_id] += G_eq
    // A[node2_id][node2_id] += G_eq
    // A[node1_id][node2_id] -= G_eq
    // A[node2_id][node1_id] -= G_eq

    if (node1_id != 0) A[node1_id][node1_id] += G_eq;
    if (node2_id != 0) A[node2_id][node2_id] += G_eq;
    if (node1_id != 0 && node2_id != 0) {
        A[node1_id][node2_id] -= G_eq;
        A[node2_id][node1_id] -= G_eq;
    }

    // Contributi al vettore lato destro B (sorgente di corrente equivalente)
    // La corrente I_srcC scorre da node1 a node2 (convenzione passiva)
    // Quindi, sottrai la corrente da node1 e aggiungi a node2
    if (node1_id != 0) B[node1_id] -= I_srcC; // Corrente che esce da node1
    if (node2_id != 0) B[node2_id] += I_srcC; // Corrente che entra in node2
}

/**
 * @brief Aggiorna lo stato interno del condensatore.
 *
 * Aggiorna la tensione e la corrente del condensatore per il prossimo passo temporale
 * utilizzando la soluzione corrente.
 *
 * @param current_solution Il vettore della soluzione corrente.
 * @param prev_solution Il vettore della soluzione precedente.
 * @param dt Il passo temporale.
 */
void Capacitor::updateState(const std::vector<double>& current_solution,
                     const std::vector<double>& prev_solution,
                     double dt) {
    if (node_ids_.size() != 2) {
        std::cerr << "Errore in updateState: Condensatore " << name_ << " si aspetta 2 ID di nodo." << std::endl;
        return;
    }
    int node1_id = node_ids_[0];
    int node2_id = node_ids_[1];

    // Controlla per indici validi
    if (node1_id < 0 || node1_id >= current_solution.size() || node2_id < 0 || node2_id >= current_solution.size()) {
        std::cerr << "Errore in updateState: Indice di nodo non valido per il Condensatore " << name_ << std::endl;
        return;
    }

    // Calcola la tensione corrente ai capi del condensatore
    double V_C_curr = (node1_id != 0 ? current_solution[node1_id] : 0.0) -
                      (node2_id != 0 ? current_solution[node2_id] : 0.0);

    // Calcola la corrente attraverso il condensatore usando la formula di integrazione trapezoidale
    // I_C(t) = (2C/dt) * (V_C(t) - V_C(t-dt)) - I_C(t-dt)
    double I_C_curr = (2.0 * capacitance / dt) * (V_C_curr - v_prev) - i_prev;

    // Aggiorna le variabili di stato per il prossimo passo temporale
    v_prev = V_C_curr;
    i_prev = I_C_curr;
}
