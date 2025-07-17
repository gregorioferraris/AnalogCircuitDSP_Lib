// components/Diode.cpp
#include "Diode.h"
#include <iostream> // Per messaggi di errore
#include <limits>   // Per std::numeric_limits (per valori minimi/massimi)

// Costante per limitare l'esponente per evitare overflow
const double MAX_EXP_ARG = 700.0; // Circa ln(DBL_MAX)

/**
 * @brief Costruttore per il componente Diode.
 * @param name Il nome univoco del diodo.
 * @param node_names_str Un vettore di stringhe contenente i nomi dei due nodi (anodo, catodo).
 * @param saturation_current La corrente di saturazione inversa (Is) in Ampere.
 * @param emission_coefficient Il coefficiente di emissione (N).
 * @param thermal_voltage La tensione termica (Vt) in Volt.
 */
Diode::Diode(const std::string& name, const std::vector<std::string>& node_names_str,
             double saturation_current, double emission_coefficient, double thermal_voltage)
    : Component(name, node_names_str),
      Is(saturation_current),
      N(emission_coefficient),
      Vt(thermal_voltage)
{
    if (node_names_str.size() != 2) {
        std::cerr << "Errore: Il diodo " << name_ << " deve avere esattamente due nodi (anodo, catodo)." << std::endl;
        // Potresti voler lanciare un'eccezione o gestire l'errore in modo più robusto.
    }
    if (Is <= 0) {
        std::cerr << "Attenzione: Corrente di saturazione (Is) non positiva per Diode " << name_ << ". Impostato a 1e-14 A." << std::endl;
        Is = 1e-14; // Valore di default ragionevole
    }
    if (N <= 0) {
        std::cerr << "Attenzione: Coefficiente di emissione (N) non positivo per Diode " << name_ << ". Impostato a 1.0." << std::endl;
        N = 1.0; // Valore di default ragionevole
    }
    if (Vt <= 0) {
        std::cerr << "Attenzione: Tensione termica (Vt) non positiva per Diode " << name_ << ". Impostato a 0.02585 V." << std::endl;
        Vt = 0.02585; // Valore di default ragionevole
    }
    std::cout << "Diode " << name_ << " inizializzato (Is: " << Is << ", N: " << N << ", Vt: " << Vt << ")." << std::endl;
}

/**
 * @brief Funzione helper per calcolare la corrente del diodo.
 * @param Vd La tensione attraverso il diodo (V_anodo - V_catodo).
 * @return La corrente che scorre attraverso il diodo.
 */
double Diode::calculateDiodeCurrent(double Vd) const {
    double arg = Vd / (N * Vt);
    // Limita l'argomento dell'esponenziale per evitare overflow
    if (arg > MAX_EXP_ARG) {
        arg = MAX_EXP_ARG;
    } else if (arg < -MAX_EXP_ARG) { // Anche per valori negativi molto grandi
        arg = -MAX_EXP_ARG;
    }
    return Is * (std::exp(arg) - 1.0);
}

/**
 * @brief Applica gli "stamps" del diodo alla matrice MNA (A) e al vettore (B).
 *
 * Questo metodo implementa il modello di linearizzazione di Newton-Raphson per il diodo.
 * Calcola la conduttanza incrementale (g) e la corrente equivalente (Ieq)
 * basandosi sulla tensione del diodo dalla soluzione corrente (x).
 *
 * @param num_total_equations Dimensione totale della matrice MNA.
 * @param dt Passo temporale (non usato direttamente per il diodo statico).
 * @param x Vettore della soluzione corrente (tensione ai nodi per il calcolo di g e Ieq).
 * @param prev_solution La soluzione dal passo temporale precedente (non usata per il diodo statico).
 * @param time Il tempo di simulazione corrente (non usato per il diodo statico).
 * @param A La matrice MNA a cui vengono applicati gli stamps.
 * @param B Il vettore lato destro MNA a cui vengono applicati gli stamps.
 */
void Diode::getStamps(
    int num_total_equations, double dt,
    const std::vector<double>& x,
    const std::vector<double>& prev_solution,
    double time,
    std::vector<std::vector<double>>& A,
    std::vector<double>& B
) {
    if (node_ids_.size() != 2) {
        std::cerr << "Errore: Diode " << name_ << " si aspetta 2 ID di nodo." << std::endl;
        return;
    }
    int anode_id = node_ids_[0];
    int cathode_id = node_ids_[1];

    // Controlla la validità degli indici
    if (anode_id < 0 || anode_id >= num_total_equations || cathode_id < 0 || cathode_id >= num_total_equations) {
        std::cerr << "Errore: Indice di nodo non valido per Diode " << name_ << std::endl;
        return;
    }

    // Ottieni la tensione del diodo dalla soluzione corrente (x)
    double V_anode = (anode_id != 0) ? x[anode_id] : 0.0;
    double V_cathode = (cathode_id != 0) ? x[cathode_id] : 0.0;
    double Vd = V_anode - V_cathode;

    // Calcola la corrente del diodo (Id) e la conduttanza incrementale (g)
    double Id = calculateDiodeCurrent(Vd);
    double g = Is / (N * Vt) * std::exp(Vd / (N * Vt)); // Derivata dId/dVd

    // Calcola la corrente equivalente per il modello di linearizzazione
    // Ieq = Id - g * Vd
    double Ieq = Id - g * Vd;

    // Applica gli stamps alla matrice A e al vettore B
    // Contributi alla matrice A (conduttanze)
    if (anode_id != 0) A[anode_id][anode_id] += g;
    if (cathode_id != 0) A[cathode_id][cathode_id] += g;
    if (anode_id != 0 && cathode_id != 0) {
        A[anode_id][cathode_id] -= g;
        A[cathode_id][anode_id] -= g;
    }

    // Contributi al vettore B (sorgenti di corrente equivalenti)
    if (anode_id != 0) B[anode_id] -= Ieq; // Corrente che esce dall'anodo
    if (cathode_id != 0) B[cathode_id] += Ieq; // Corrente che entra nel catodo
}

/**
 * @brief Aggiorna lo stato interno del diodo.
 *
 * Per un diodo ideale, non c'è uno stato interno che evolve nel tempo.
 * Questo metodo è lasciato vuoto per questo modello.
 *
 * @param current_solution Il vettore della soluzione corrente (non usato).
 * @param prev_solution Il vettore della soluzione precedente (non usato).
 * @param dt Il passo temporale (non usato).
 */
void Diode::updateState(const std::vector<double>& current_solution,
                     const std::vector<double>& prev_solution,
                     double dt) {
    // Nessuno stato interno da aggiornare per un diodo ideale.
    // Questo metodo è intenzionalmente lasciato vuoto.
}
