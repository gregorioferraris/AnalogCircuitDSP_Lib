// components/Attenuator.h
#ifndef ATTENUATOR_H
#define ATTENUATOR_H

#include "Component.h" // Include la classe base Component
#include <Eigen/Dense> // Per le operazioni con matrici e vettori Eigen
#include <string>      // Per std::string

/**
 * @class Attenuator
 * @brief Rappresenta un attenuatore ideale in una simulazione di circuito.
 *
 * Questa classe modella un attenuatore ideale come un dispositivo a quattro terminali.
 * La tensione di uscita è una versione scalata della tensione di ingresso, definita da un
 * fattore di attenuazione. Questo modello si comporta come un generatore di tensione
 * controllato in tensione (VCVS) con guadagno costante.
 *
 * Si assume che l'attenuatore abbia impedenza di ingresso infinita (non assorbe corrente
 * dai nodi di ingresso) e impedenza di uscita nulla (può fornire qualsiasi corrente
 * necessaria per mantenere la tensione di uscita).
 *
 * Introduce una variabile ausiliaria nella matrice MNA per la corrente che scorre
 * attraverso il lato di uscita dell'attenuatore.
 *
 * L'equazione implementata è:
 * V_output = attenuation_factor * V_input
 * dove V_input = V(node1_in) - V(node2_in)
 * e V_output = V(node1_out) - V(node2_out)
 */
class Attenuator : public Component {
public:
    // Nomi dei nodi per il lato di uscita
    std::string node1_out; ///< Nome del primo nodo sul lato di uscita (terminale positivo).
    std::string node2_out; ///< Nome del secondo nodo sul lato di uscita (terminale negativo).

    double attenuation_factor; ///< Il fattore di attenuazione (K) dell'attenuatore (V_out = K * V_in).

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
     * @param tolerance_percent Percentuale di tolleranza opzionale (non usata direttamente per l'attenuatore ideale).
     */
    Attenuator(const std::string& name,
               const std::string& node1_in, const std::string& node2_in,
               const std::string& node1_out, const std::string& node2_out,
               double attenuation_factor, double tolerance_percent = 0.0);

    /**
     * @brief Crea una copia profonda dell'oggetto Attenuator.
     * @return Un puntatore a un nuovo oggetto Attenuator, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new Attenuator(*this); }

    /**
     * @brief Applica gli "stamps" dell'attenuatore ideale alla matrice MNA (A) e al vettore (b).
     *
     * Questo metodo introduce una variabile ausiliaria per la corrente di uscita.
     * Gli stamps impongono la relazione di tensione V_output = attenuation_factor * V_input.
     *
     * @param A La matrice MNA a cui vengono applicati gli stamps.
     * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
     * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo (non usata per questo componente statico).
     * @param prev_solution La soluzione dal passo temporale precedente (non usata per questo componente statico).
     * @param time Il tempo di simulazione corrente (non usata per questo componente statico).
     * @param dt La dimensione del passo temporale (non usata per questo componente statico).
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Aggiorna lo stato interno dell'attenuatore.
     *
     * Essendo un attenuatore ideale un componente statico e senza perdite (il suo comportamento
     * dipende solo dalle tensioni e correnti istantanee, non dagli stati passati),
     * questo metodo non esegue alcun aggiornamento di stato. È fornito per soddisfare
     * l'interfaccia della classe base Component.
     *
     * @param v_curr La tensione corrente attraverso il componente (non usata).
     * @param i_curr La corrente corrente che attraversa il componente (non usata).
     */
    void updateState(double v_curr, double i_curr) override;
};

#endif // ATTENUATOR_H
