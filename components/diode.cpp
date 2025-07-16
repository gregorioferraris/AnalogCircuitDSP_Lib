#include "Diode.h"
#include <iostream>

// Costruttore
Diode::Diode(const std::string& name, const std::string& anode_node, const std::string& cathode_node,
             double Is, double N, double Vt)
    : Component(name, anode_node, cathode_node), Is_(Is), N_(N), Vt_(Vt)
{
    // Il costruttore della classe base Component gestisce l'inizializzazione dei nodi.
}

// Implementazione di calculateCurrent
double Diode::calculateCurrent(double Vd) const {
    double arg = Vd / (N_ * Vt_);
    
    // Limiti numerici per prevenire overflow/underflow in std::exp
    // Equivalentemente a np.exp in Python
    const double EXP_MAX_ARG = 700.0; // Valore approssimativo per std::exp che non causi overflow per double
    const double EXP_MIN_ARG = -70.0; // Valore approssimativo per std::exp che non causi underflow significativo

    if (arg > EXP_MAX_ARG) {
        return Is_ * (std::exp(EXP_MAX_ARG) - 1.0);
    } else if (arg < EXP_MIN_ARG) {
        // In inversa profonda, la corrente tende a -Is
        // Usiamo un piccolo valore per stabilità, ma è molto vicino a -Is
        return -Is_; // Potrebbe anche essere Is_ * (std::exp(arg) - 1.0) se arg è solo < -70, ma non così piccolo da essere 0
                     // Per arg < -70, exp(arg) è estremamente vicino a 0, quindi -Is è una buona approssimazione.
    } else {
        return Is_ * (std::exp(arg) - 1.0);
    }
}

// Implementazione di calculateConductance
double Diode::calculateConductance(double Vd) const {
    double arg = Vd / (N_ * Vt_);

    const double EXP_MAX_ARG = 700.0;
    const double SMALL_CONDUCTANCE = 1e-9; // Piccola conduttanza per stabilità in inversa

    if (arg > EXP_MAX_ARG) {
        // In forte conduzione, la conduttanza è molto alta
        return Is_ / (N_ * Vt_) * std::exp(EXP_MAX_ARG);
    } else if (arg < EXP_MIN_ARG) { // Utilizziamo lo stesso limite di EXP_MIN_ARG
        // In inversa profonda, la conduttanza è molto piccola (quasi zero)
        return SMALL_CONDUCTANCE; // Non esattamente zero per evitare divisioni per zero o problemi numerici
    } else {
        return Is_ / (N_ * Vt_) * std::exp(arg);
    }
}

// Implementazione di getStamps per il diodo
void Diode::getStamps(Eigen::MatrixXd& stamp_A, Eigen::VectorXd& stamp_B,
                   int num_total_equations, double dt,
                   const Eigen::VectorXd& current_solution_guess,
                   const Eigen::VectorXd& prev_solution,
                   double time)
{
    // Verifica e ridimensiona/azzerra matrici se necessario (solitamente fatto dal solutore)
    if (stamp_A.rows() != num_total_equations || stamp_A.cols() != num_total_equations) {
        stamp_A = Eigen::MatrixXd::Zero(num_total_equations, num_total_equations);
    }
    if (stamp_B.size() != num_total_equations) {
        stamp_B = Eigen::VectorXd::Zero(num_total_equations);
    }

    // Ottieni gli ID numerici dei nodi
    if (node_ids.size() != 2 || node_ids[0] == -1 || node_ids[1] == -1) {
        std::cerr << "Errore: ID dei nodi non validi per il diodo " << name << std::endl;
        return;
    }

    int anode_id = node_ids[0];
    int cathode_id = node_ids[1];

    // Ottieni le tensioni correnti dai nodi della soluzione "current_solution_guess"
    // current_solution_guess contiene le tensioni nodali e possibilmente correnti aggiuntive.
    // Assumiamo che i primi 'num_nodes' elementi di current_solution_guess siano le tensioni nodali.
    double V_anode = (anode_id != 0) ? current_solution_guess(anode_id) : 0.0;
    double V_cathode = (cathode_id != 0) ? current_solution_guess(cathode_id) : 0.0;

    double Vd = V_anode - V_cathode; // Tensione ai capi del diodo
    double Id = calculateCurrent(Vd);     // Corrente del diodo
    double Gd = calculateConductance(Vd); // Conduttanza dinamica del diodo

    // Applica lo "stamp" del diodo al sistema di equazioni Newton-Raphson.
    // Il contributo è dato dalla conduttanza dinamica e dalla corrente.
    // L'equazione per un diodo in MNA è:
    // I_diode(Vd) = Gd * Vd - I_diode(Vd) + Gd * Vd
    // Questa forma è equivalente a:
    // I_diode(Vd) = Gd * Vd + (Id - Gd * Vd)
    // Dove (Id - Gd * Vd) è il termine di corrente equivalente per il metodo Newton-Raphson.

    // Contributi alla matrice Jacobiana (stamp_A):
    // La conduttanza Gd contribuisce ai termini diagonali e incrociati.
    // Se Anodo e Catodo non sono massa (0)
    if (anode_id != 0) {
        stamp_A(anode_id, anode_id) += Gd;
    }
    if (cathode_id != 0) {
        stamp_A(cathode_id, cathode_id) += Gd;
    }
    if (anode_id != 0 && cathode_id != 0) {
        stamp_A(anode_id, cathode_id) -= Gd;
        stamp_A(cathode_id, anode_id) -= Gd;
    }

    // Contributi al vettore del lato destro (stamp_B) (corrente di mancato equilibrio):
    // Il termine Id - Gd * Vd viene sottratto dalle equazioni nodali.
    // Ricorda che la corrente esce dal nodo positivo (anodo) ed entra nel negativo (catodo).
    double current_contribution = Id - Gd * Vd;

    if (anode_id != 0) {
        stamp_B(anode_id) -= current_contribution;
    }
    if (cathode_id != 0) {
        stamp_B(cathode_id) += current_contribution;
    }
}
