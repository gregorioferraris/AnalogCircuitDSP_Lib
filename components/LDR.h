// components/LDR.h
#ifndef LDR_H
#define LDR_H

#include "Component.h" // Include la classe base Component
#include <Eigen/Dense> // Per le operazioni con matrici e vettori Eigen
#include <string>      // Per std::string
#include <cmath>       // Per std::pow

/**
 * @class LDR
 * @brief Rappresenta un Resistore Fotoelettrico (Light Dependent Resistor).
 *
 * La resistenza dell'LDR varia in base all'intensità luminosa secondo la formula:
 * R_LDR = R0 * (L / L0)^(-gamma)
 *
 * Dove:
 * - R_LDR: Resistenza attuale dell'LDR.
 * - R0: Resistenza di riferimento (a L0 lux).
 * - L: Intensità luminosa attuale (lux).
 * - L0: Intensità luminosa di riferimento (lux).
 * - gamma: Coefficiente di sensibilità alla luce.
 */
class LDR : public Component {
public:
    // Parametri del modello LDR
    double R0;       ///< Resistenza di riferimento (Ohm) a L0 lux.
    double L0;       ///< Intensità luminosa di riferimento (lux).
    double gamma;    ///< Coefficiente di sensibilità alla luce (adimensionale).

    // Variabile di stato per l'intensità luminosa attuale
    double current_light_intensity; ///< Intensità luminosa attuale (lux).

    /**
     * @brief Costruttore per il componente LDR.
     *
     * Inizializza un LDR con il suo nome, i due nodi collegati e i parametri del modello.
     *
     * @param name Il nome univoco dell'LDR.
     * @param node1 Il nome del primo nodo.
     * @param node2 Il nome del secondo nodo.
     * @param R0 La resistenza di riferimento (Ohm) a L0 lux. Deve essere > 0.
     * @param L0 L'intensità luminosa di riferimento (lux). Deve essere > 0.
     * @param gamma Il coefficiente di sensibilità alla luce. Deve essere > 0.
     * @param initial_light_intensity L'intensità luminosa iniziale (lux). Deve essere >= 0.
     * @param tolerance_percent Percentuale di tolleranza opzionale (non usata direttamente per questo modello).
     */
    LDR(const std::string& name,
        const std::string& node1, const std::string& node2,
        double R0, double L0, double gamma, double initial_light_intensity = 10.0,
        double tolerance_percent = 0.0);

    /**
     * @brief Crea una copia profonda dell'oggetto LDR.
     * @return Un puntatore a un nuovo oggetto LDR, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new LDR(*this); }

    /**
     * @brief Applica gli "stamps" del modello LDR alla matrice MNA (A) e al vettore (b).
     *
     * Calcola la resistenza attuale dell'LDR in base a `current_light_intensity`
     * e applica lo stamp di un resistore.
     *
     * @param A La matrice MNA a cui vengono applicati gli stamps.
     * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
     * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo (non usata direttamente).
     * @param prev_solution La soluzione dal passo temporale precedente (non usata direttamente).
     * @param time Il tempo di simulazione corrente (non usato direttamente).
     * @param dt La dimensione del passo temporale (non usata direttamente).
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Imposta l'intensità luminosa attuale dell'LDR.
     *
     * Questo metodo permette di aggiornare l'intensità luminosa esternamente,
     * influenzando la resistenza dell'LDR nel passo di simulazione successivo.
     *
     * @param light_intensity La nuova intensità luminosa in lux. Deve essere >= 0.
     */
    void setLightIntensity(double light_intensity);

    /**
     * @brief Aggiorna lo stato interno del componente (non usato per l'LDR).
     *
     * Per l'LDR, la resistenza è determinata dall'intensità luminosa esterna,
     * non dalle sue tensioni o correnti precedenti. Questo metodo è vuoto.
     *
     * @param v_curr La tensione corrente attraverso il componente (non usata).
     * @param i_curr La corrente che attraversa il componente (non usata).
     */
    void updateState(double v_curr, double i_curr) override {
        // LDR resistance is determined by external light intensity, not its own voltage/current.
        // No internal state update based on v_curr or i_curr is needed for this model.
    }
};

#endif // LDR_H
