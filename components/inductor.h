// components/Inductor.h
#ifndef INDUCTOR_H
#define INDUCTOR_H

#include "Component.h"
#include <Eigen/Dense>

class Inductor : public Component {
public:
    double L;           // Valore effettivo dell'induttanza
    double i_prev;      // Corrente attraverso l'induttore al passo precedente
    double v_prev;      // Tensione ai capi dell'induttore al passo precedente

    Inductor(const std::string& name, const std::string& node1, const std::string& node2,
             double inductance_nominal, double tolerance_percent = 0.0);
    
    Component* clone() const override { return new Inductor(*this); } // Implementazione del clone

    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;
    
    void updateState(double v_curr, double i_curr) override;
};

#endif // INDUCTOR_H
