// components/DynamicMicrophone.cpp
#include "DynamicMicrophone.h"
#include <iostream> // Per potenziali messaggi di debug

/**
 * @brief Costruttore per il componente DynamicMicrophone.
 * @param name Il nome univoco del microfono.
 * @param node1 Il nome del nodo del terminale positivo.
 * @param node2 Il nome del nodo del terminale negativo.
 * @param v_func Una funzione che restituisce la tensione di uscita del microfono a un dato tempo.
 * @param tolerance_percent Percentuale di tolleranza opzionale (non usata direttamente per la sorgente ideale).
 */
DynamicMicrophone::DynamicMicrophone(const std::string& name, const std::string& node1, const std::string& node2,
                                     std::function<double(double)> v_func, double tolerance_percent)
    : Component(name, node1, node2, tolerance_percent), voltage_function(std::move(v_func))
{
    // Un microfono dinamico, modellato come una sorgente di tensione ideale, richiede una variabile ausiliaria
    // per la sua corrente di ramo nella formulazione MNA.
    // La classe base Component dovrebbe avere un meccanismo per indicare questa necessità,
    // e il simulatore assegnerà l'indice effettivo.
    // Assumiamo che il costruttore della classe base Component o una chiamata successiva
    // gestisca l'indicizzazione della variabile ausiliaria.
    setRequiresAuxiliaryVariable(true);
}

/**
 * @brief Crea una copia profonda dell'oggetto DynamicMicrophone.
 * @return Un puntatore a un nuovo oggetto DynamicMicrophone, che è una copia dell'istanza corrente.
 */
Component* DynamicMicrophone::clone() const {
    // Durante la clonazione, assicurarsi che la std::function sia copiata correttamente.
    // std::function ha un costruttore di copia.
    return new DynamicMicrophone(*this);
}

/**
 * @brief Applica gli "stamps" del microfono dinamico (come sorgente di tensione)
 * alla matrice MNA (A) e al vettore (b).
 *
 * Questo componente introduce una variabile ausiliaria per la sua corrente di ramo.
 * Gli stamps impongono il vincolo di tensione $V(node1) - V(node2) = \text{voltage\_function(time)}$.
 *
 * @param A La matrice MNA a cui vengono applicati gli stamps.
 * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
 * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo (non usata direttamente qui).
 * @param prev_solution La soluzione dal passo temporale precedente (non usata qui).
 * @param time Il tempo di simulazione corrente.
 * @param dt La dimensione del passo temporale (non usata qui).
 */
void DynamicMicrophone::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Ottieni gli indici globali per i due nodi collegati.
    int idx1 = getNodeIndex(node1); // Terminale positivo
    int idx2 = getNodeIndex(node2); // Terminale negativo

    // Ottieni l'indice per la variabile ausiliaria di corrente di questa sorgente di tensione.
    // Questo indice deve essere stato assegnato dal simulatore.
    int aux_idx = getAuxiliaryVariableIndex();

    // Controlla se gli indici sono validi. Se uno è -1, significa che un nodo o una variabile ausiliaria
    // non è stata correttamente registrata/trovata nella mappatura globale.
    if (idx1 == -1 || idx2 == -1 || aux_idx == -1) {
        std::cerr << "Errore: Indice di nodo o ausiliario non valido per DynamicMicrophone " << name << std::endl;
        // A seconda della strategia di gestione degli errori, potrebbe lanciare un'eccezione o ritornare.
        return;
    }

    // Ottieni la tensione di uscita corrente dalla funzione del microfono.
    double V_mic = voltage_function(time);

    // Applica gli stamps per una sorgente di tensione indipendente:
    // 1. Conservazione della corrente al nodo1: $+I_{mic}$ (corrente che esce dal nodo1 attraverso il microfono)
    //    Questo significa che la variabile di corrente $I_{mic}$ contribuisce alla somma delle correnti al nodo1.
    //    In MNA, questo si traduce in $A(\text{idx1}, \text{aux\_idx}) += 1$.
    A(idx1, aux_idx) += 1.0;

    // 2. Conservazione della corrente al nodo2: $-I_{mic}$ (corrente che entra nel nodo2 attraverso il microfono)
    //    In MNA, questo si traduce in $A(\text{idx2}, \text{aux\_idx}) -= 1$.
    A(idx2, aux_idx) -= 1.0;

    // 3. Vincolo di tensione: $V(node1) - V(node2) = V_{mic}$
    //    Questo vincolo è posto nella riga corrispondente alla variabile ausiliaria di corrente.
    //    $A(\text{aux\_idx}, \text{idx1}) += 1$ (coefficiente per $V(node1)$)
    //    $A(\text{aux\_idx}, \text{idx2}) -= 1$ (coefficiente per $V(node2)$)
    //    $b(\text{aux\_idx}) = V_{mic}$ (la tensione della sorgente)
    A(aux_idx, idx1) += 1.0;
    A(aux_idx, idx2) -= 1.0;
    b(aux_idx) += V_mic; // Aggiungi la tensione della sorgente al RHS
}

/**
 * @brief Aggiorna lo stato interno del microfono dinamico.
 *
 * Essendo una sorgente di tensione ideale, l'uscita del microfono è determinata unicamente
 * dalla funzione di tensione fornita in funzione del tempo e non ha uno stato interno
 * che evolve in base alla propria tensione o corrente.
 * Questo metodo è intenzionalmente lasciato vuoto.
 *
 * @param v_curr La tensione corrente attraverso il componente (non usata).
 * @param i_curr La corrente corrente che attraversa il componente (non usata).
 */
void DynamicMicrophone::updateState(double v_curr, double i_curr) {
    // Nessuno stato interno da aggiornare per una sorgente di tensione ideale.
}
