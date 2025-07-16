// components/RibbonTweeter.h
#ifndef RIBBON_TWEETER_H
#define RIBBON_TWEETER_H

#include "Component.h" // Include la classe base Component
#include <string>

/**
 * @brief Rappresenta un componente Ribbon Tweeter in una simulazione di circuito.
 *
 * Questo componente modella un tweeter a nastro come un carico elettrico,
 * costituito da una resistenza in serie e un'induttanza.
 *
 * Ha due terminali elettrici: un nodo di ingresso/positivo e un nodo di uscita/negativo.
 * Per l'analisi transitoria, l'induttore è modellato utilizzando il suo modello di accompagnamento.
 * Questo componente introduce una variabile ausiliaria per la corrente che lo attraversa.
 */
class RibbonTweeter : public Component {
public:
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
    RibbonTweeter(const std::string& name,
                  const std::string& positive_node,
                  const std::string& negative_node,
                  double ribbon_resistance,
                  double ribbon_inductance,
                  double tolerance_percent = 0.0);

    /**
     * @brief Applica gli "stamps" del RibbonTweeter alla matrice MNA (A) e al vettore (b).
     *
     * Questo metodo implementa il modello di accompagnamento per la combinazione in serie di una
     * resistenza e un induttore, rappresentando il carico del tweeter.
     * Aggiunge una variabile ausiliaria per rappresentare la corrente che scorre attraverso il tweeter.
     *
     * @param A La matrice MNA a cui vengono applicati gli "stamps".
     * @param b Il vettore lato destro MNA a cui vengono applicati gli "stamps".
     * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti dei rami.
     * @param prev_solution La soluzione del passo temporale precedente, utilizzata per il modello di accompagnamento dell'induttore.
     * @param time Il tempo di simulazione corrente (non direttamente usato per questo componente passivo, ma richiesto dall'interfaccia).
     * @param dt La dimensione del passo temporale.
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

private:
    double R_ribbon; // Resistenza del nastro (Ohm)
    double L_ribbon; // Induttanza del nastro (Henry)

    // Nomi dei nodi per chiarezza (già memorizzati nella classe base Component, ma utili come copie locali)
    std::string pos_node_name;
    std::string neg_node_name;
};

#endif // RIBBON_TWEETER_H
