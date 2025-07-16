// components/Capacitor.cpp
#include "Capacitor.h"
// #include "../utils/Random.h" // Per applyTolerance - Assumiamo che non sia più necessario o gestito altrove
#include <stdexcept>
#include <iostream> // Per i messaggi di errore

/**
 * @brief Costruttore per il componente Capacitor.
 *
 * @param name Il nome univoco del condensatore.
 * @param node_names_str Un vettore di stringhe contenente i nomi dei due nodi collegati.
 * @param capacitance_nominal La capacità nominale in Farad.
 */
Capacitor::Capacitor(const std::string& name, const std::vector<std::string>& node_names_str,
                     double capacitance_nominal)
    : Component(name, node_names_str), v_prev_(0.0), i_prev_(0.0)
{
    // Non usiamo applyTolerance qui per semplicità, se non è definito o non necessario.
    // capacitance = applyTolerance(capacitance_nominal, tolerance_percent);
    capacitance_ = capacitance_nominal;

    if (capacitance_ <= 0) {
        std::cerr << "Attenzione: Capacità non positiva per Capacitor " << name_ << ". Impostato a 1e-12F." << std::endl;
        capacitance_ = 1e-12; // Valore di default piccolo ma positivo
    }
    std::cout << "Capacitor " << name_ << " inizializzato con capacità: " << capacitance_ << " F." << std::endl;
}

/**
 * @brief Applica gli "stamps" del condensatore alla matrice MNA (A) e al vettore (B).
 *
 * Questo metodo utilizza il modello del compagno di Eulero all'indietro per il condensatore
 * per applicare i suoi contributi in base allo stato del passo temporale precedente.
 *
 * @param num_total_equations Dimensione totale della matrice MNA.
 * @param dt Passo temporale. Cruciale per i modelli del compagno.
 * @param x Vettore della soluzione corrente (non usato per lo stamping, ma per updateState).
 * @param prev_solution Vettore della soluzione al passo temporale precedente (usato per i modelli del compagno).
 * @param time Tempo attuale della simulazione.
 * @param A Riferimento alla matrice MNA.
 * @param B Riferimento al vettore delle sorgenti (RHS).
 */
void Capacitor::getStamps(
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
        std::cerr << "Errore: Capacitor " << name_ << " si aspetta 2 ID di nodo, ma ne ha " << node_ids_.size() << std::endl;
        return;
    }
    int node1_id = node_ids_[0];
    int node2_id = node_ids_[1];

    // Controlla la validità degli indici
    if (node1_id < 0 || node1_id >= num_total_equations || node2_id < 0 || node2_id >= num_total_equations) {
        std::cerr << "Errore: Indice di nodo non valido per Capacitor " << name_ << std::endl;
        return;
    }

    if (dt <= 0) {
        std::cerr << "Errore: Il passo temporale (dt) deve essere positivo per lo stamping del condensatore " << name_ << "." << std::endl;
        return;
    }

    // G_eq = C / dt (per Eulero in avanti) o 2*C / dt (per Trapezoidale o Eulero all'indietro con corrente media)
    // Usiamo il modello del compagno di Eulero all'indietro per coerenza con altri componenti dinamici.
    // I_C(t) = C * dV_C/dt  =>  I_C(t) = C * (V_C(t) - V_C(t-dt)) / dt
    // Questo è equivalente a una conduttanza G_eq = C / dt in parallelo con una sorgente di corrente I_src = (C / dt) * V_C(t-dt)
    double G_eq = capacitance_ / dt;

    // Contributi alla matrice MNA (parte dipendente dalla tensione attuale)
    // Gli stamps sono come quelli di un resistore con conduttanza G_eq.
    if (node1_id != 0) A[node1_id][node1_id] += G_eq;
    if (node2_id != 0) A[node2_id][node2_id] += G_eq;
    if (node1_id != 0 && node2_id != 0) {
        A[node1_id][node2_id] -= G_eq;
        A[node2_id][node1_id] -= G_eq;
    }

    // Contributi al vettore RHS (parte dipendente dallo stato precedente)
    // V_C_prev è la tensione ai capi del condensatore al passo temporale precedente.
    // Il nodo 0 è il ground.
    double V_C_prev = ((node1_id == 0) ? 0.0 : prev_solution[node1_id]) -
                      ((node2_id == 0) ? 0.0 : prev_solution[node2_id]);
    double I_src = G_eq * V_C_prev;

    // La corrente I_src fluisce da node2 a node1.
    // Quindi, aggiungiamo I_src a B[node1_id] e sottraiamo da B[node2_id].
    if (node1_id != 0) B[node1_id] += I_src;
    if (node2_id != 0) B[node2_id] -= I_src;
}

/**
 * @brief Aggiorna le variabili di stato interne (v_prev_, i_prev_) per il prossimo passo temporale.
 *
 * Questo metodo viene chiamato dopo che il sistema MNA è stato risolto per il passo temporale corrente.
 * Aggiorna la tensione ai capi del condensatore (`v_prev_`) per il prossimo passo.
 * La corrente (`i_prev_`) può essere calcolata se necessaria, ma per il modello di Eulero all'indietro
 * per i condensatori, solo la tensione precedente è strettamente necessaria per lo stamp.
 *
 * @param current_solution Il vettore della soluzione corrente (contiene le tensioni ai nodi).
 * @param prev_solution Il vettore della soluzione precedente (non usato direttamente qui per l'aggiornamento).
 * @param dt Il passo temporale.
 */
void Capacitor::updateState(const std::vector<double>& current_solution,
                         const std::vector<double>& prev_solution,
                         double dt) {
    // Ottieni gli indici globali per i due nodi collegati.
    int node1_id = node_ids_[0];
    int node2_id = node_ids_[1];

    // Calcola la tensione corrente ai capi del condensatore.
    // Il nodo 0 è il ground.
    v_prev_ = ((node1_id == 0) ? 0.0 : current_solution[node1_id]) -
              ((node2_id == 0) ? 0.0 : current_solution[node2_id]);

    // La corrente i_prev_ può essere calcolata se necessaria per altri scopi,
    // ma non è strettamente richiesta per lo stamp del condensatore con Eulero all'indietro.
    // i_prev_ = capacitance_ * (v_prev_ - (prev_solution[node1_id] - prev_solution[node2_id])) / dt;
}
