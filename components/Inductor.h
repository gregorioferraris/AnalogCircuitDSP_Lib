// components/Inductor.h
#ifndef INDUCTOR_H
#define INDUCTOR_H

#include "Component.h"
#include <Eigen/Dense>

class Inductor : public Component {
public:
    double L;           // Effective inductance value
    double i_prev;      // Current through the inductor at the previous step
    double v_prev;      // Voltage across the inductor at the previous step

    Inductor(const std::string& name, const std::string& node1, const std::string& node2,
             double inductance_nominal, double tolerance_percent = 0.0);
    
    Component* clone() const override { return new Inductor(*this); } // Clone implementation

    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;
    
    void updateState(double v_curr, double i_curr) override;
};

#endif // INDUCTOR_H
