// components/Pentode.h
#ifndef PENTODE_H
#define PENTODE_H

#include "Component.h" // Include la classe base Component
#include <Eigen/Dense>   // Per le operazioni con matrici e vettori Eigen
#include <string>        // Per std::string
#include <vector>        // Per std::vector
#include <map>           // Per std::map (per getNodeIndex)

/**
 * @class Pentode
 * @brief Rappresenta un componente pentodo per simulazioni di circuiti MNA.
 *
 * Questa classe modella un pentodo con i suoi cinque terminali principali:
 * Griglia di Controllo (Control Grid), Placca (Plate), Catodo (Cathode),
 * Griglia Schermo (Screen Grid), e Griglia Soppressore (Suppressor Grid).
 * Implementa i metodi per applicare gli "stamps" alla matrice MNA e al vettore RHS,
 * utilizzando un modello non lineare semplificato adatto all'analisi transitoria
 * con iterazione di Newton-Raphson.
 */
class Pentode : public Component {
public:
    // Nomi dei nodi per i terminali del pentodo
    std::string node_control_grid;
    std::string node_plate;
    std::string node_cathode;
    std::string node_screen_grid;
    std::string node_suppressor_grid;

    // Parametri del modello semplificato del pentodo
    double K_p_;        ///< Parametro di transconduttanza (simile a Kp in MOSFET)
    double V_cutoff_;   ///< Tensione di cutoff della griglia di controllo
    double lambda_;     ///< Parametro di modulazione della lunghezza del canale (per dipendenza dalla placca)
    double mu_sg_;      ///< Fattore di amplificazione della griglia schermo (quanto Vsg influenza Iplate rispetto a Vcg)
    double I_s_ratio_;  ///< Rapporto tra corrente di griglia schermo e corrente di placca (semplificato)

    // Variabili di stato per l'iterazione di Newton-Raphson (tensioni del passo precedente)
    // Queste sono le tensioni dei nodi calcolate all'iterazione precedente del Newton-Raphson.
    double prev_Vcg_;   ///< Tensione precedente della griglia di controllo
    double prev_Vp_;    ///< Tensione precedente della placca
    double prev_Vsg_;   ///< Tensione precedente della griglia schermo
    double prev_Vsup_;  ///< Tensione precedente della griglia soppressore

    /**
     * @brief Costruttore per il componente Pentode.
     * @param name Il nome univoco del pentodo.
     * @param control_grid_node Il nome del nodo della griglia di controllo.
     * @param plate_node Il nome del nodo della placca.
     * @param cathode_node Il nome del nodo del catodo.
     * @param screen_grid_node Il nome del nodo della griglia schermo.
     * @param suppressor_grid_node Il nome del nodo della griglia soppressore.
     * @param K_p Parametro di transconduttanza.
     * @param V_cutoff Tensione di cutoff.
     * @param lambda Parametro lambda per la modulazione della lunghezza del canale.
     * @param mu_sg Fattore di amplificazione della griglia schermo.
     * @param I_s_ratio Rapporto corrente griglia schermo/placca.
     * @param tolerance_percent Percentuale di tolleranza opzionale.
     */
    Pentode(const std::string& name,
            const std::string& control_grid_node, const std::string& plate_node,
            const std::string& cathode_node, const std::string& screen_grid_node,
            const std::string& suppressor_grid_node,
            double K_p, double V_cutoff, double lambda, double mu_sg, double I_s_ratio,
            double tolerance_percent = 0.0);

    /**
     * @brief Crea una copia profonda dell'oggetto Pentode.
     * @return Un puntatore a un nuovo oggetto Pentode, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new Pentode(*this); }

    /**
     * @brief Applica gli "stamps" del pentodo alla matrice MNA (A) e al vettore (b).
     *
     * Questo metodo implementa un modello non lineare semplificato per il pentodo,
     * calcolando le correnti di placca e griglia schermo e le loro conduttanze dinamiche
     * per l'iterazione di Newton-Raphson.
     *
     * @param A La matrice MNA a cui vengono applicati gli stamps.
     * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
     * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo.
     * @param prev_solution La soluzione dal passo temporale precedente (non usata direttamente per il modello statico non lineare).
     * @param time Il tempo di simulazione corrente.
     * @param dt La dimensione del passo temporale.
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Aggiorna lo stato interno del pentodo.
     *
     * Questo metodo memorizza le tensioni dei nodi attuali come "precedenti"
     * per la prossima iterazione di Newton-Raphson nel metodo getStamps.
     *
     * @param v_curr La tensione corrente attraverso il componente (non usata direttamente).
     * @param i_curr La corrente corrente che attraversa il componente (non usata direttamente).
     */
    void updateState(double v_curr, double i_curr) override;
};

#endif // PENTODE_H
