// components/Capacitor.h
#ifndef CAPACITOR_H
#define CAPACITOR_H

#include "Component.h"
#include <Eigen/Dense>

class Capacitor : public Component {
public:
    double capacitance; // Valore effettivo della capacità
    double v_prev;      // Tensione ai capi del condensatore al passo precedente
    double i_prev;      // Corrente attraverso il condensatore al passo precedente

    Capacitor(const std::string& name, const std::string& node1, const std::string& node2,
              double capacitance_nominal, double tolerance_percent = 0.0);
    
    Component* clone() const override { return new Capacitor(*this); } // Implementazione del clone

    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;
    
    void updateState(double v_curr, double i_curr) override;
};

#endif // CAPACITOR_H
