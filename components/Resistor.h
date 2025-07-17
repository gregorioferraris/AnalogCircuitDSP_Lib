// components/Resistor.h
#ifndef RESISTOR_H
#define RESISTOR_H

#include "Component.h" // Include la classe base Component
#include <string>        // Per std::string
#include <vector>        // Per std::vector<std::string>

/**
 * @class Resistor
 * @brief Rappresenta un componente resistore in una simulazione di circuito.
 *
 * Questa classe modella un resistore lineare con una resistenza fissa.
 * Implementa i metodi per applicare i "stamps" alla matrice MNA e al vettore RHS.
 */
class Resistor : public Component {
public:
    double resistance; // Valore della resistenza in Ohm

    /**
     * @brief Costruttore per il componente Resistor.
     * @param name Il nome univoco del resistore.
     * @param node_names_str Un vettore di stringhe contenente i nomi dei due nodi collegati.
     * @param resistance_nominal Il valore nominale della resistenza in Ohm.
     */
    Resistor(const std::string& name, const std::vector<std::string>& node_names_str,
             double resistance_nominal);

    /**
     * @brief Crea una copia profonda dell'oggetto Resistor.
     * @return Un puntatore a un nuovo oggetto Resistor, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new Resistor(*this); }

    /**
     * @brief Applica gli "stamps" del resistore alla matrice MNA (A) e al vettore (B).
     *
     * Per un resistore, gli stamps sono costanti e lineari, basati sulla legge di Ohm.
     *
     * @param num_total_equations Dimensione totale della matrice MNA.
     * @param dt Passo temporale (non usato per i resistori statici).
     * @param x Vettore della soluzione corrente (non usato per i resistori lineari).
     * @param prev_solution La soluzione dal passo temporale precedente (non usata per i resistori statici).
     * @param time Il tempo di simulazione corrente (non usato per i resistori statici).
     * @param A La matrice MNA a cui vengono applicati gli stamps.
     * @param B Il vettore lato destro MNA a cui vengono applicati gli stamps.
     */
    void getStamps(
        int num_total_equations, double dt,
        const std::vector<double>& x,
        const std::vector<double>& prev_solution,
        double time,
        std::vector<std::vector<double>>& A,
        std::vector<double>& B
    ) override;

    /**
     * @brief Aggiorna lo stato interno del resistore.
     *
     * Per un resistore ideale, non c'è uno stato interno che evolve nel tempo.
     * Questo metodo è lasciato vuoto per questo modello.
     *
     * @param current_solution Il vettore della soluzione corrente (non usato).
     * @param prev_solution Il vettore della soluzione precedente (non usato).
     * @param dt Il passo temporale (non usato).
     */
    void updateState(const std::vector<double>& current_solution,
                     const std::vector<double>& prev_solution,
                     double dt) override;
};

#endif // RESISTOR_H
