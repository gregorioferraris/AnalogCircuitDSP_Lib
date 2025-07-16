// components/RibbonTweeter.cpp
#include "RibbonTweeter.h"
#include <iostream> // Per output di debug e avvisi

/**
 * @brief Costruttore per il componente RibbonTweeter.
 *
 * Inizializza un tweeter a nastro con il suo nome, due nodi di connessione,
 * resistenza e induttanza del nastro.
 *
 * @param name Il nome univoco del tweeter.
 * @param positive_node Il nome del nodo di ingresso/positivo.
 * @param negative_node Il nome del nodo di uscita/negativo.
 * @param ribbon_resistance La resistenza in serie del nastro (Ohm).
 * @param ribbon_inductance L'induttanza in serie del nastro (Henry).
 * @param tolerance_percent Tolleranza opzionale in percentuale (ereditata da Component).
 */
RibbonTweeter::RibbonTweeter(const std::string& name,
                             const std::string& positive_node,
                             const std::string& negative_node,
                             double ribbon_resistance,
                             double ribbon_inductance,
                             double tolerance_percent)
    : Component(name, positive_node, negative_node, tolerance_percent), // Passa i due nodi principali alla classe base Component
      R_ribbon(ribbon_resistance), L_ribbon(ribbon_inductance),
      pos_node_name(positive_node), neg_node_name(negative_node)
{
    // Un tweeter a nastro, modellato come una resistenza e un induttore in serie,
    // introduce una variabile ausiliaria per la sua corrente di ramo in MNA.
    setNumAuxiliaryVariables(1);

    // Validazione di base per i parametri
    if (R_ribbon < 0 || L_ribbon < 0) {
        std::cerr << "Attenzione: RibbonTweeter " << name << " ha resistenza (" << R_ribbon << " Ohm) o induttanza (" << L_ribbon << " Henry) negative. Questo potrebbe portare a simulazioni instabili." << std::endl;
    }
    std::cout << "RibbonTweeter " << name << " inizializzato. Connesso tra " << positive_node << " e " << negative_node << "." << std::endl;
}

/**
 * @brief Applica gli "stamps" del RibbonTweeter alla matrice MNA (A) e al vettore (b).
 *
 * Questo metodo implementa il modello di accompagnamento per la combinazione in serie di una
 * resistenza e un induttore.
 *
 * L'induttore L è modellato utilizzando il modello di accompagnamento di Eulero all'indietro:
 * V_L(t) = (L/dt) * I_L(t) - (L/dt) * I_L(t-dt)
 *
 * L'equazione totale del ramo (dal nodo positivo al nodo negativo, con corrente I_tweeter):
 * V_pos - V_neg = I_tweeter(t) * R_ribbon + V_L(t)
 * Sostituisci V_L(t):
 * V_pos - V_neg = I_tweeter(t) * R_ribbon + (L_ribbon/dt) * I_tweeter(t) - (L_ribbon/dt) * I_tweeter(t-dt)
 *
 * Riorganizza nella forma MNA (termini con incognite a sinistra, termini noti a destra):
 * V_pos - V_neg - (R_ribbon + L_ribbon/dt) * I_tweeter(t) = -(L_ribbon/dt) * I_tweeter(t-dt)
 *
 * Sia R_eff = R_ribbon + L_ribbon/dt
 * Sia V_eff_source = -(L_ribbon/dt) * I_tweeter(t-dt)
 *
 * Gli "stamps" MNA sono:
 * 1. KCL al nodo positivo: -I_tweeter
 * 2. KCL al nodo negativo: +I_tweeter
 * 3. Equazione costitutiva (riga ausiliaria): V_pos - V_neg - R_eff * I_tweeter = V_eff_source
 *
 * @param A La matrice MNA a cui vengono applicati gli "stamps".
 * @param b Il vettore lato destro MNA a cui vengono applicati gli "stamps".
 * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti dei rami (non direttamente usata per gli "stamps" di questo modello lineare).
 * @param prev_solution La soluzione del passo temporale precedente, usata per il termine storico del modello di accompagnamento dell'induttore.
 * @param time Il tempo di simulazione corrente.
 * @param dt La dimensione del passo temporale.
 */
void RibbonTweeter::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Ottieni gli indici globali (base 1) per i nodi positivo e negativo.
    // Assumendo che getNodeIndex(node_name) restituisca l'indice base 1 (0 per la massa).
    int idx_p_1based = getNodeIndex(pos_node_name);
    int idx_n_1based = getNodeIndex(neg_node_name);

    // Ottieni l'indice globale (base 1) per la variabile di corrente ausiliaria introdotta da questo componente.
    // Assumendo che getAuxiliaryVariableStartIndex() restituisca l'indice base 1 della prima (e unica)
    // variabile ausiliaria per questo componente tweeter.
    int idx_tweeter_current_1based = getAuxiliaryVariableStartIndex();

    // --- Calcola i termini per il modello di accompagnamento ---

    // Calcola la resistenza efficace per il modello di accompagnamento dell'induttore (Eulero all'indietro)
    // Per l'analisi DC (dt=0), un induttore agisce come un cortocircuito, quindi R_L_eff = 0.
    // Per l'analisi transitoria (dt > 0), R_L_eff = L_ribbon / dt.
    double R_L_eff = 0.0;
    if (dt > 1e-12) { // Controlla se dt è effettivamente non zero per l'analisi transitoria
        R_L_eff = L_ribbon / dt;
    }
    // Se dt è 0, R_L_eff rimane 0, il che modella correttamente l'induttore come un cortocircuito in DC.

    // Ottieni la corrente precedente attraverso il ramo del tweeter da prev_solution.
    // Questo è il termine storico I_tweeter(t-dt) per il modello di accompagnamento dell'induttore.
    // Se prev_solution non è ancora popolato (es. primo passo temporale), assumi che la corrente precedente sia 0.
    double I_tweeter_prev = 0.0;
    if (prev_solution.size() > idx_tweeter_current_1based) {
        I_tweeter_prev = prev_solution(idx_tweeter_current_1based);
    }

    // Calcola la resistenza totale efficace in serie del ramo del tweeter
    // Questo include la resistenza del nastro e la resistenza efficace dal modello di accompagnamento dell'induttore.
    double R_eff = R_ribbon + R_L_eff;

    // Calcola la sorgente di tensione efficace per il modello di accompagnamento.
    // Questo termine rappresenta il contributo del termine storico dell'induttore.
    // V_eff_source = -(L_ribbon/dt) * I_tweeter(t-dt)
    double V_eff_source = -R_L_eff * I_tweeter_prev; // Il termine di sorgente dal modello di accompagnamento dell'induttore

    // --- Applica gli "stamps" alla matrice MNA A e al vettore b ---

    // 1. Equazione KCL per il nodo positivo (idx_p_1based):
    // La corrente I_tweeter esce dal nodo positivo. Quindi, il suo coefficiente è -1.
    // A(riga_idx, col_idx) += valore
    A(idx_p_1based, idx_tweeter_current_1based) += -1.0;

    // 2. Equazione KCL per il nodo negativo (idx_n_1based):
    // La corrente I_tweeter entra nel nodo negativo. Quindi, il suo coefficiente è +1.
    A(idx_n_1based, idx_tweeter_current_1based) += 1.0;

    // 3. Equazione costitutiva per il ramo del tweeter (riga ausiliaria idx_tweeter_current_1based):
    // L'equazione è: V_pos - V_neg - R_eff * I_tweeter = V_eff_source
    // Coefficiente per V_pos è +1
    A(idx_tweeter_current_1based, idx_p_1based) += 1.0;
    // Coefficiente per V_neg è -1
    A(idx_tweeter_current_1based, idx_n_1based) += -1.0;
    // Coefficiente per I_tweeter è -R_eff
    A(idx_tweeter_current_1based, idx_tweeter_current_1based) += -R_eff;

    // Aggiungi la sorgente di tensione efficace al vettore lato destro (RHS) b.
    // Questo termine è un valore noto per il passo temporale corrente.
    b(idx_tweeter_current_1based) += V_eff_source;
}
