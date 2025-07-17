// components/Diode.h
#ifndef DIODE_H
#define DIODE_H

#include "Component.h" // Include la classe base Component
#include <string>        // Per std::string
#include <vector>        // Per std::vector
#include <cmath>         // Per std::exp

/**
 * @class Diode
 * @brief Rappresenta un componente diodo in una simulazione di circuito.
 *
 * Questa classe modella un diodo ideale con una corrente di saturazione e un coefficiente di emissione.
 * È un componente non lineare e richiede un solutore iterativo (es. Newton-Raphson)
 * per una simulazione accurata.
 */
class Diode : public Component {
public:
    double Is; // Corrente di saturazione (A)
    double N;  // Coefficiente di emissione (o fattore di idealità)
    double Vt; // Tensione termica (V) - tipicamente 0.02585V a 300K

    /**
     * @brief Costruttore per il componente Diode.
     * @param name Il nome univoco del diodo.
     * @param node_names_str Un vettore di stringhe contenente i nomi dei due nodi (anodo, catodo).
     * @param saturation_current La corrente di saturazione inversa (Is) in Ampere.
     * @param emission_coefficient Il coefficiente di emissione (N).
     * @param thermal_voltage La tensione termica (Vt) in Volt.
     */
    Diode(const std::string& name, const std::vector<std::string>& node_names_str,
          double saturation_current, double emission_coefficient, double thermal_voltage = 0.02585);

    /**
     * @brief Crea una copia profonda dell'oggetto Diode.
     * @return Un puntatore a un nuovo oggetto Diode, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new Diode(*this); }

    /**
     * @brief Applica gli "stamps" del diodo alla matrice MNA (A) e al vettore (B).
     *
     * Per un diodo, gli stamps sono non lineari e dipendono dalla tensione corrente.
     * Questo metodo calcola la conduttanza incrementale (g) e la corrente equivalente (Ieq)
     * per il modello di linearizzazione di Newton-Raphson.
     *
     * @param num_total_equations Dimensione totale della matrice MNA.
     * @param dt Passo temporale (non usato direttamente per il diodo statico).
     * @param x Vettore della soluzione corrente (tensione ai nodi per il calcolo di g e Ieq).
     * @param prev_solution La soluzione dal passo temporale precedente (non usata per il diodo statico).
     * @param time Il tempo di simulazione corrente (non usato per il diodo statico).
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
     * @brief Aggiorna lo stato interno del diodo.
     *
     * Per un diodo ideale, non c'è uno stato interno che evolve nel tempo.
     * Questo metodo è lasciato vuoto per questo modello.
     *
     * @param current_solution Il vettore della soluzione corrente (non usato).
     * @param prev_solution Il vettore della soluzione precedente (non usato).
     * @param dt Il passo temporale (non usato).
     */
    void updateState(const std::vector<double>& current_solution,
                     const std::vector<double>& prev_solution,
                     double dt) override;

private:
    // Funzione helper per calcolare la corrente del diodo
    double calculateDiodeCurrent(double Vd) const;
};

#endif // DIODE_H
