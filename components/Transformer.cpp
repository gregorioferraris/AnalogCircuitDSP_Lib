// components/Transformer.cpp
#include "Transformer.h"
#include <iostream> // Per i messaggi di errore

/**
 * @brief Costruttore per il componente Transformer.
 *
 * Inizializza un trasformatore ideale con il suo nome, i quattro nodi collegati
 * e il suo rapporto di spire.
 *
 * @param name Il nome univoco del trasformatore.
 * @param node_names_str Un vettore di stringhe contenente i nomi dei nodi nell'ordine:
 * [node1_p, node2_p, node1_s, node2_s] (primario positivo, primario negativo, secondario positivo, secondario negativo).
 * @param turns_ratio Il rapporto di spire (Np/Ns) del trasformatore. Deve essere diverso da zero.
 */
Transformer::Transformer(const std::string& name,
                         const std::vector<std::string>& node_names_str,
                         double turns_ratio)
    : Component(name, node_names_str), // Passa tutti i nomi dei nodi al costruttore della classe base
      turns_ratio_(turns_ratio)
{
    // Un trasformatore ideale richiede due variabili ausiliarie per le sue correnti primaria (Ip) e secondaria (Is).
    // Questo è gestito dal MnaSolver che assegna gli indici delle variabili ausiliarie.
    // Non è necessario chiamare setNumAuxiliaryVariables() qui, ma il MnaSolver deve sapere
    // quanti ne servono per ogni componente.
    // Per ora, assumiamo che MnaSolver gestisca questo basandosi sul tipo di componente
    // o che un metodo getNumAuxiliaryVariables() sia implementato nella classe base Component.
    if (turns_ratio_ == 0.0) {
        std::cerr << "Attenzione: Transformer " << name_ << " ha un rapporto di spire pari a zero. Questo può portare a divisioni per zero." << std::endl;
    }
    std::cout << "Transformer " << name_ << " inizializzato con rapporto di spire: " << turns_ratio_ << std::endl;
}

/**
 * @brief Applica gli "stamps" del trasformatore ideale alla matrice MNA (A) e al vettore (B).
 *
 * Questo metodo introduce due variabili ausiliarie per le correnti primaria (Ip) e secondaria (Is).
 * Applica le relazioni di tensione e corrente di un trasformatore ideale.
 *
 * @param num_total_equations Dimensione totale della matrice MNA.
 * @param dt Passo temporale per la simulazione transitoria (non usato per questo componente statico).
 * @param x Vettore della soluzione corrente (per ottenere le tensioni dei nodi).
 * @param prev_solution Vettore della soluzione al passo temporale precedente (non usato).
 * @param time Tempo attuale della simulazione (non usato).
 * @param A Riferimento alla matrice MNA.
 * @param B Riferimento al vettore delle sorgenti (RHS).
 */
void Transformer::getStamps(
    int num_total_equations, double dt,
    const std::vector<double>& x,
    const std::vector<double>& prev_solution,
    double time,
    std::vector<std::vector<double>>& A,
    std::vector<double>& B
) {
    // Ottieni gli indici globali per tutti e quattro i nodi collegati.
    // Assumiamo che node_ids_ contenga [node1_p, node2_p, node1_s, node2_s]
    // nell'ordine: Primario Positivo (0), Primario Negativo (1), Secondario Positivo (2), Secondario Negativo (3)
    if (node_ids_.size() != 4) {
        std::cerr << "Errore: Transformer " << name_ << " si aspetta 4 ID di nodo, ma ne ha " << node_ids_.size() << std::endl;
        return;
    }

    int idx_p1 = node_ids_[0];
    int idx_p2 = node_ids_[1];
    int idx_s1 = node_ids_[2];
    int idx_s2 = node_ids_[3];

    // Ottieni gli indici per le variabili ausiliarie (correnti primaria e secondaria).
    // Questi indici devono essere gestiti dal MnaSolver.
    // Come per l'Attenuatore, useremo un indice fittizio temporaneo.
    // L'MnaSolver dovrebbe assegnare questi indici in modo consecutivo.
    // Assumiamo che le variabili ausiliarie siano aggiunte dopo tutti i nodi.
    // Per un trasformatore, ci sono 2 variabili ausiliarie.
    // L'indice della prima variabile ausiliaria per questo componente sarà `num_total_equations - 2`
    // e la seconda `num_total_equations - 1` se è l'ultimo componente che aggiunge variabili ausiliarie.
    // Questo è un placeholder e deve essere gestito correttamente dal MnaSolver.
    // Idealmente, Component dovrebbe avere un metodo `getAuxiliaryVariableStartIndex()`.
    // Per ora, useremo una convenzione che dovrà essere allineata con l'MnaSolver.
    // Se il MnaSolver assegna `component_id_` come indice di partenza per le variabili ausiliarie,
    // allora `aux_idx_Ip = component_id_` e `aux_idx_Is = component_id_ + 1`.
    // Questo richiede che MnaSolver sia consapevole di quanti aux_vars ogni componente ha.
    // Per ora, userò un indice basato su una stima del numero totale di variabili ausiliarie.
    // Questo è un *workaround* temporaneo.
    int aux_idx_Ip = num_total_equations - 2; // Placeholder per la corrente primaria
    int aux_idx_Is = num_total_equations - 1; // Placeholder per la corrente secondaria

    // Assicurati che turns_ratio non sia zero per evitare divisioni per zero.
    if (turns_ratio_ == 0.0) {
        std::cerr << "Errore: Transformer " << name_ << " ha un rapporto di spire pari a zero, impossibile applicare gli stamps." << std::endl;
        return;
    }

    // --- Applica gli stamps per le equazioni del trasformatore ideale ---

    // Equazione 1: V_primario - turns_ratio * V_secondario = 0
    // (V(idx_p1) - V(idx_p2)) - turns_ratio * (V(idx_s1) - V(idx_s2)) = 0
    // Questa equazione di vincolo va nella riga corrispondente alla variabile ausiliaria della corrente primaria (Ip).
    A[aux_idx_Ip][idx_p1] += 1.0;
    A[aux_idx_Ip][idx_p2] -= 1.0;
    A[aux_idx_Ip][idx_s1] -= turns_ratio_;
    A[aux_idx_Ip][idx_s2] += turns_ratio_;
    // B[aux_idx_Ip] rimane 0.0 poiché è un'equazione omogenea.

    // Equazione 2: Ip + (1/turns_ratio) * Is = 0
    // Questa equazione va nella riga corrispondente alla variabile ausiliaria della corrente secondaria (Is).
    A[aux_idx_Is][aux_idx_Ip] += 1.0;
    A[aux_idx_Is][aux_idx_Is] += 1.0 / turns_ratio_;
    // B[aux_idx_Is] rimane 0.0 poiché è un'equazione omogenea.

    // --- Applica i contributi di corrente alle equazioni dei nodi ---

    // Corrente lato primario: Ip entra in idx_p1 ed esce da idx_p2
    A[idx_p1][aux_idx_Ip] += 1.0;
    A[idx_p2][aux_idx_Ip] -= 1.0;

    // Corrente lato secondario: Is entra in idx_s1 ed esce da idx_s2
    A[idx_s1][aux_idx_Is] += 1.0;
    A[idx_s2][aux_idx_Is] -= 1.0;
}

/**
 * @brief Aggiorna lo stato interno del trasformatore.
 *
 * Essendo un trasformatore ideale un componente statico e senza perdite (il suo comportamento
 * dipende solo dalle tensioni e correnti istantanee, non dagli stati passati),
 * questo metodo non esegue alcun aggiornamento di stato. È fornito per soddisfare
 * l'interfaccia della classe base Component, garantendo che tutti i metodi virtuali siano implementati.
 *
 * @param current_solution Il vettore della soluzione corrente (non usato).
 * @param prev_solution Il vettore della soluzione precedente (non usato).
 * @param dt Il passo temporale (non usato).
 */
void Transformer::updateState(const std::vector<double>& current_solution,
                         const std::vector<double>& prev_solution,
                         double dt) {
    // Nessuno stato interno da aggiornare per un modello di trasformatore ideale.
    // Questo metodo è intenzionalmente lasciato vuoto.
}
