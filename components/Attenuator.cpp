// components/Attenuator.cpp
#include "Attenuator.h"
#include <iostream> // Per messaggi di errore

/**
 * @brief Costruttore per il componente Attenuator.
 *
 * Inizializza un attenuatore ideale con il suo nome, i quattro nodi collegati
 * e il suo fattore di attenuazione.
 *
 * @param name Il nome univoco dell'attenuatore.
 * @param node1_in Il nome del primo nodo sul lato di ingresso (terminale positivo).
 * @param node2_in Il nome del secondo nodo sul lato di ingresso (terminale negativo).
 * @param node1_out Il nome del primo nodo sul lato di uscita (terminale positivo).
 * @param node2_out Il nome del secondo nodo sul lato di uscita (terminale negativo).
 * @param attenuation_factor Il fattore di attenuazione (es. 0.5 per 6dB di attenuazione).
 * @param tolerance_percent Percentuale di tolleranza opzionale.
 */
Attenuator::Attenuator(const std::string& name,
                       const std::string& node1_in, const std::string& node2_in,
                       const std::string& node1_out, const std::string& node2_out,
                       double attenuation_factor, double tolerance_percent)
    // Passa i nodi di ingresso al costruttore della classe base Component.
    // node1 e node2 della classe base sono i nodi di ingresso.
    : Component(name, node1_in, node2_in, tolerance_percent),
      node1_out(node1_out), node2_out(node2_out), attenuation_factor(attenuation_factor)
{
    // Un attenuatore ideale (modellato come un VCVS) richiede una variabile ausiliaria
    // per la sua corrente di uscita.
    setNumAuxiliaryVariables(1); // Richiede una variabile ausiliaria
}

/**
 * @brief Applica gli "stamps" dell'attenuatore ideale alla matrice MNA (A) e al vettore (b).
 *
 * Questo metodo introduce una variabile ausiliaria per la corrente di uscita.
 * Gli stamps impongono la relazione di tensione V_output = attenuation_factor * V_input.
 *
 * @param A La matrice MNA a cui vengono applicati gli stamps.
 * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
 * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo (non usata).
 * @param prev_solution La soluzione dal passo temporale precedente (non usata).
 * @param time Il tempo di simulazione corrente (non usata).
 * @param dt La dimensione del passo temporale (non usata).
 */
void Attenuator::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Ottieni gli indici globali per tutti e quattro i nodi.
    // Si assume che getNodeIndex sia un metodo della classe base Component
    // (o accessibile tramite Circuit) che mappa correttamente i nomi dei nodi
    // ai loro indici corrispondenti nella matrice MNA.
    int idx_in1 = getNodeIndex(node1);  // Nodo di ingresso positivo (Component::node1)
    int idx_in2 = getNodeIndex(node2);  // Nodo di ingresso negativo (Component::node2)
    int idx_out1 = getNodeIndex(node1_out); // Nodo di uscita positivo
    int idx_out2 = getNodeIndex(node2_out); // Nodo di uscita negativo

    // Ottieni l'indice per la variabile ausiliaria (corrente di uscita).
    int aux_idx = getAuxiliaryVariableStartIndex(); // Ottiene l'indice di partenza per le variabili ausiliarie

    // Controlla se gli indici sono validi. Se un indice è -1, significa che un nodo
    // o una variabile ausiliaria non è stata correttamente registrata/trovata nella mappatura globale.
    if (idx_in1 == -1 || idx_in2 == -1 || idx_out1 == -1 || idx_out2 == -1 || aux_idx == -1) {
        std::cerr << "Errore: Indice di nodo o ausiliario non valido per Attenuator " << name << std::endl;
        return;
    }

    // --- Applica gli stamps per l'equazione dell'attenuatore ideale (VCVS) ---
    // Equazione: V_output - attenuation_factor * V_input = 0
    // (V(node1_out) - V(node2_out)) - attenuation_factor * (V(node1_in) - V(node2_in)) = 0

    // Questa equazione di vincolo va nella riga corrispondente alla variabile ausiliaria.
    A(aux_idx, idx_out1) += 1.0;
    A(aux_idx, idx_out2) -= 1.0;
    A(aux_idx, idx_in1) -= attenuation_factor;
    A(aux_idx, idx_in2) += attenuation_factor;
    // b(aux_idx) rimane 0.0 poiché è un'equazione omogenea.

    // --- Applica i contributi di corrente alle equazioni dei nodi (lato di uscita) ---
    // La variabile ausiliaria rappresenta la corrente che esce da node1_out ed entra in node2_out.
    // Questo è il comportamento di una sorgente di tensione controllata.
    A(idx_out1, aux_idx) += 1.0;
    A(idx_out2, aux_idx) -= 1.0;

    // Nota: Un VCVS ideale ha impedenza di ingresso infinita, quindi non vengono applicati
    // stamps di corrente ai nodi di ingresso (idx_in1, idx_in2) da parte di questo componente.
}

/**
 * @brief Aggiorna lo stato interno dell'attenuatore.
 *
 * Essendo un attenuatore ideale (modellato come un VCVS statico), non ha uno stato interno
 * da aggiornare. Questo metodo è intenzionalmente lasciato vuoto per soddisfare l'interfaccia
 * della classe base Component.
 *
 * @param v_curr La tensione corrente attraverso il componente (non usata).
 * @param i_curr La corrente corrente che attraversa il componente (non usata).
 */
void Attenuator::updateState(double v_curr, double i_curr) {
    // Nessuno stato interno da aggiornare per un attenuatore ideale.
}
