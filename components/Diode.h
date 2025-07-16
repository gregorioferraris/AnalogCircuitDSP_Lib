#ifndef DIODE_H
#define DIODE_H

#include <string>
#include <vector>
#include <cmath> // Per std::exp
#include <limits> // Per std::numeric_limits

#include <Eigen/Dense> // Per Eigen::MatrixXd e Eigen::VectorXd

// Assicurati che Component sia già definito (es. in Component.h)
// Nel nostro caso, lo abbiamo già fatto in CurrentSource.h, quindi potresti
// includere direttamente quel file o creare un Component.h separato.
// Per ora, lo includiamo da CurrentSource.h per avere la base Component.
#include "CurrentSource.h" // Include la definizione della classe Component

// --- Classe Diode ---
class Diode : public Component {
public:
    // Costruttore
    Diode(const std::string& name, const std::string& anode_node, const std::string& cathode_node,
          double Is = 1e-14, double N = 1.0, double Vt = 0.0258);

    // Implementazione del metodo virtuale getStamps dalla classe base Component
    // Per i componenti non lineari, questo metodo contribuisce alla matrice Jacobiana (stamp_A)
    // e al vettore di "mancato equilibrio" (stamp_B) per Newton-Raphson.
    void getStamps(Eigen::MatrixXd& stamp_A, Eigen::VectorXd& stamp_B,
                   int num_total_equations, double dt,
                   const Eigen::VectorXd& current_solution_guess,
                   const Eigen::VectorXd& prev_solution,
                   double time) override;

    // Calcola la corrente attraverso il diodo data la tensione ai suoi capi (Vd = V_anode - V_cathode)
    double calculateCurrent(double Vd) const;

    // Calcola la conduttanza dinamica del diodo (dI/dVd)
    double calculateConductance(double Vd) const;

private:
    double Is_; // Corrente di saturazione inversa
    double N_;  // Fattore di idealità
    double Vt_; // Tensione termica (kT/q)
};

#endif // DIODE_H
