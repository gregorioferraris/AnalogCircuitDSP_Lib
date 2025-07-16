// components/RibbonMicrophone.cpp
#include "RibbonMicrophone.h"
#include <iostream> // Per output di debug e avvisi

/**
 * @brief Costruttore per il componente RibbonMicrophone.
 *
 * Inizializza un microfono a nastro con il suo nome, due nodi di uscita,
 * resistenza di uscita, induttanza di uscita e una funzione che fornisce
 * la pressione sonora istantanea in funzione del tempo.
 *
 * @param name Il nome univoco del microfono.
 * @param positive_node Il nome del nodo di uscita positivo.
 * @param negative_node Il nome del nodo di uscita negativo.
 * @param output_resistance La resistenza di uscita in serie del microfono (Ohm).
 * @param output_inductance L'induttanza di uscita in serie del microfono (Henry).
 * @param sensitivity La sensibilità del microfono (Volt/Pascal).
 * @param sound_pressure_function Una std::function che prende il tempo corrente (double)
 * e restituisce la pressione sonora istantanea (Pascal).
 * @param tolerance_percent Tolleranza opzionale in percentuale (ereditata da Component).
 */
RibbonMicrophone::RibbonMicrophone(const std::string& name,
                                   const std::string& positive_node,
                                   const std::string& negative_node,
                                   double output_resistance,
                                   double output_inductance,
                                   double sensitivity,
                                   std::function<double(double time)> sound_pressure_function,
                                   double tolerance_percent)
    : Component(name, positive_node, negative_node, tolerance_percent), // Passa i due nodi principali alla classe base Component
      R_out(output_resistance), L_out(output_inductance),
      sensitivity(sensitivity),
      V_sound_pressure_func(sound_pressure_function),
      pos_node_name(positive_node), neg_node_name(negative_node)
{
    // Un microfono a nastro, modellato come una sorgente di tensione con impedenza in serie,
    // introduce una variabile ausiliaria per la sua corrente di ramo in MNA.
    setNumAuxiliaryVariables(1);

    // Validazione di base per i parametri
    if (R_out < 0 || L_out < 0 || sensitivity < 0) {
        std::cerr << "Attenzione: RibbonMicrophone " << name << " ha resistenza di uscita (" << R_out << " Ohm), induttanza (" << L_out << " Henry) o sensibilità (" << sensitivity << " V/Pa) negative. Questo potrebbe portare a simulazioni instabili." << std::endl;
    }
    std::cout << "RibbonMicrophone " << name << " inizializzato. Uscita tra " << positive_node << " e " << negative_node << "." << std::endl;
}

/**
 * @brief Applica gli "stamps" del RibbonMicrophone alla matrice MNA (A) e al vettore (b).
 *
 * Questo metodo implementa il modello di accompagnamento per la combinazione in serie di una
 * sorgente di tensione variabile nel tempo, una resistenza e un induttore.
 *
 * L'induttore L è modellato utilizzando il modello di accompagnamento di Eulero all'indietro:
 * V_L(t) = (L/dt) * I_L(t) - (L/dt) * I_L(t-dt)
 *
 * L'equazione totale del ramo (dal nodo positivo al nodo negativo, con corrente I_mic):
 * V_pos - V_neg = V_mic_generated(t) + I_mic(t) * R_out + V_L(t)
 * Sostituisci V_L(t):
 * V_pos - V_neg = V_mic_generated(t) + I_mic(t) * R_out + (L_out/dt) * I_mic(t) - (L_out/dt) * I_mic(t-dt)
 *
 * Riorganizza nella forma MNA (termini con incognite a sinistra, termini noti a destra):
 * V_pos - V_neg - (R_out + L_out/dt) * I_mic(t) = V_mic_generated(t) - (L_out/dt) * I_mic(t-dt)
 *
 * Sia R_total_eff = R_out + L_out/dt
 * Sia V_eff_source = V_mic_generated(t) - (L_out/dt) * I_mic(t-dt)
 *
 * Gli "stamps" MNA sono:
 * 1. KCL al nodo positivo: -I_mic
 * 2. KCL al nodo negativo: +I_mic
 * 3. Equazione costitutiva (riga ausiliaria): V_pos - V_neg - R_total_eff * I_mic = V_eff_source
 *
 * @param A La matrice MNA a cui vengono applicati gli "stamps".
 * @param b Il vettore lato destro MNA a cui vengono applicati gli "stamps".
 * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti dei rami (non direttamente usata per gli "stamps" di questo modello lineare).
 * @param prev_solution La soluzione del passo temporale precedente, usata per il termine storico del modello di accompagnamento dell'induttore.
 * @param time Il tempo di simulazione corrente.
 * @param dt La dimensione del passo temporale.
 */
void RibbonMicrophone::getStamps(
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
    // variabile ausiliaria per questo componente microfono.
    int idx_mic_current_1based = getAuxiliaryVariableStartIndex();

    // --- Calcola i termini per il modello di accompagnamento ---

    // Calcola la pressione sonora istantanea
    double sound_pressure_at_t = V_sound_pressure_func(time);
    // Calcola la tensione generata dal microfono
    double V_mic_generated = sensitivity * sound_pressure_at_t;

    // Calcola la resistenza efficace per il modello di accompagnamento dell'induttore (Eulero all'indietro)
    // Per l'analisi DC (dt=0), un induttore agisce come un cortocircuito, quindi R_L_eff = 0.
    // Per l'analisi transitoria (dt > 0), R_L_eff = L_out / dt.
    double R_L_eff = 0.0;
    if (dt > 1e-12) { // Controlla se dt è effettivamente non zero per l'analisi transitoria
        R_L_eff = L_out / dt;
    }
    // Se dt è 0, R_L_eff rimane 0, il che modella correttamente l'induttore come un cortocircuito in DC.

    // Ottieni la corrente precedente attraverso il ramo del microfono da prev_solution.
    // Questo è il termine storico I_mic(t-dt) per il modello di accompagnamento dell'induttore.
    // Se prev_solution non è ancora popolato (es. primo passo temporale), assumi che la corrente precedente sia 0.
    double I_mic_prev = 0.0;
    if (prev_solution.size() > idx_mic_current_1based) {
        I_mic_prev = prev_solution(idx_mic_current_1based);
    }

    // Calcola la resistenza totale efficace in serie del ramo del microfono
    // Questo include la resistenza di uscita e la resistenza efficace dal modello di accompagnamento dell'induttore.
    double R_total_eff = R_out + R_L_eff;

    // Calcola la sorgente di tensione efficace per il modello di accompagnamento.
    // Questo include la tensione generata dal suono e il termine storico dall'induttore.
    // V_eff_source = V_mic_generated(t) - (L_out/dt) * I_mic(t-dt)
    double V_L_companion_source = R_L_eff * I_mic_prev; // La parte di sorgente di tensione del modello di accompagnamento dell'induttore
    double V_eff_source = V_mic_generated - V_L_companion_source;

    // --- Applica gli "stamps" alla matrice MNA A e al vettore b ---

    // 1. Equazione KCL per il nodo positivo (idx_p_1based):
    // La corrente I_mic esce dal nodo positivo. Quindi, il suo coefficiente è -1.
    // A(riga_idx, col_idx) += valore
    A(idx_p_1based, idx_mic_current_1based) += -1.0;

    // 2. Equazione KCL per il nodo negativo (idx_n_1based):
    // La corrente I_mic entra nel nodo negativo. Quindi, il suo coefficiente è +1.
    A(idx_n_1based, idx_mic_current_1based) += 1.0;

    // 3. Equazione costitutiva per il ramo del microfono (riga ausiliaria idx_mic_current_1based):
    // L'equazione è: V_pos - V_neg - R_total_eff * I_mic = V_eff_source
    // Coefficiente per V_pos è +1
    A(idx_mic_current_1based, idx_p_1based) += 1.0;
    // Coefficiente per V_neg è -1
    A(idx_mic_current_1based, idx_n_1based) += -1.0;
    // Coefficiente per I_mic è -R_total_eff
    A(idx_mic_current_1based, idx_mic_current_1based) += -R_total_eff;

    // Aggiungi la sorgente di tensione efficace al vettore lato destro (RHS) b.
    // Questo termine è un valore noto per il passo temporale corrente.
    b(idx_mic_current_1based) += V_eff_source;
}
