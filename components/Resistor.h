// components/Resistor.h
#ifndef RESISTOR_H
#define RESISTOR_H

#include "Component.h" // Include la classe base Component
#include <Eigen/Dense>   // Per le operazioni con matrici e vettori Eigen
#include <string>        // Per std::string

/**
 * @class Resistor
 * @brief Rappresenta un componente resistore in una simulazione di circuito.
 *
 * Questa classe modella un resistore lineare con una resistenza fissa.
 * Implementa i metodi per applicare i "stamps" alla matrice MNA e al vettore RHS.
 */
class Resistor : public Component {
public:
    double R; // Valore della resistenza in Ohm

    /**
     * @brief Costruttore per il componente Resistor.
     * @param name Il nome univoco del resistore.
     * @param node1 Il nome del primo nodo (terminale).
     * @param node2 Il nome del secondo nodo (terminale).
     * @param resistance Il valore della resistenza in Ohm.
     * @param tolerance_percent Percentuale di tolleranza opzionale per il valore del componente.
     */
    Resistor(const std::string& name, const std::string& node1, const std::string& node2,
             double resistance, double tolerance_percent = 0.0);

    /**
     * @brief Crea una copia profonda dell'oggetto Resistor.
     * @return Un puntatore a un nuovo oggetto Resistor, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new Resistor(*this); }

    /**
     * @brief Applica gli "stamps" del resistore alla matrice MNA (A) e al vettore (b).
     *
     * Per un resistore, gli stamps sono costanti e lineari, basati sulla legge di Ohm.
     *
     * @param A La matrice MNA a cui vengono applicati gli stamps.
     * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
     * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo (non usata per i resistori lineari).
     * @param prev_solution La soluzione dal passo temporale precedente (non usata per i resistori statici).
     * @param time Il tempo di simulazione corrente (non usato per i resistori statici).
     * @param dt La dimensione del passo temporale (non usata per i resistori statici).
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Aggiorna lo stato interno del resistore.
     *
     * Per un resistore ideale, non c'è uno stato interno che evolve nel tempo.
     * Questo metodo è lasciato vuoto per questo modello.
     *
     * @param v_curr La tensione corrente attraverso il componente.
     * @param i_curr La corrente corrente che attraversa il componente.
     */
    void updateState(double v_curr, double i_curr) override;
};

#endif // RESISTOR_H
