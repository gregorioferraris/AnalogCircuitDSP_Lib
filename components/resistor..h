// components/Resistor.h
#ifndef RESISTOR_H
#define RESISTOR_H

#include "Component.h"
#include <Eigen/Dense> // Per Eigen::MatrixXd, Eigen::VectorXd

class Resistor : public Component {
public:
    double resistance; // Valore effettivo della resistenza

    Resistor(const std::string& name, const std::string& node1, const std::string& node2,
             double resistance_nominal, double tolerance_percent = 0.0);
    
    Component* clone() const override { return new Resistor(*this); } // Implementazione del clone

    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;
};

#endif // RESISTOR_H
