// components/Capacitor.h
#ifndef CAPACITOR_H
#define CAPACITOR_H

#include "Component.h" // Include la classe base Component
#include <string>        // Per std::string
#include <vector>        // Per std::vector

/**
 * @class Capacitor
 * @brief Rappresenta un componente condensatore in una simulazione di circuito.
 *
 * Questa classe modella un condensatore lineare. Implementa i metodi per applicare
 * i "stamps" alla matrice MNA e al vettore RHS usando il metodo di integrazione
 * trapezoidale per l'analisi transitoria.
 */
class Capacitor : public Component {
public:
    double capacitance; // Valore della capacità in Farad

    // Variabili di stato per l'integrazione trapezoidale
    double v_prev; // Tensione ai capi del condensatore al passo temporale precedente
    double i_prev; // Corrente attraverso il condensatore al passo temporale precedente

    /**
     * @brief Costruttore per il componente Capacitor.
     * @param name Il nome univoco del condensatore.
     * @param node_names_str Un vettore di stringhe contenente i nomi dei due nodi collegati.
     * @param capacitance_nominal Il valore nominale della capacità in Farad.
     */
    Capacitor(const std::string& name, const std::vector<std::string>& node_names_str,
              double capacitance_nominal);

    /**
     * @brief Crea una copia profonda dell'oggetto Capacitor.
     * @return Un puntatore a un nuovo oggetto Capacitor, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new Capacitor(*this); }

    /**
     * @brief Applica gli "stamps" del condensatore alla matrice MNA (A) e al vettore (B).
     *
     * Utilizza il metodo di integrazione trapezoidale per modellare il condensatore
     * in analisi transitoria.
     *
     * @param num_total_equations Dimensione totale della matrice MNA.
     * @param dt Passo temporale.
     * @param x Vettore della soluzione corrente (tensione ai nodi e correnti di ramo).
     * @param prev_solution La soluzione dal passo temporale precedente.
     * @param time Il tempo di simulazione corrente.
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
     * @brief Aggiorna lo stato interno del condensatore.
     *
     * Aggiorna la tensione e la corrente del condensatore per il passo temporale successivo.
     *
     * @param current_solution Il vettore della soluzione corrente.
     * @param prev_solution Il vettore della soluzione precedente.
     * @param dt Il passo temporale.
     */
    void updateState(const std::vector<double>& current_solution,
                     const std::vector<double>& prev_solution,
                     double dt) override;
};

#endif // CAPACITOR_H
