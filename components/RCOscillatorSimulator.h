// File: RCOscillatorSimulator.h
#ifndef RC_OSCILLATOR_SIMULATOR_H
#define RC_OSCILLATOR_SIMULATOR_H

/**
 * @class RCOscillatorSimulator
 * @brief Simula un semplice oscillatore a rilassamento RC.
 *
 * Questa classe modella il comportamento di un condensatore che si carica e si scarica
 * attraverso un resistore, controllato da un comparatore con isteresi.
 * L'output è la tensione del condensatore e lo stato logico del comparatore.
 *
 * Utilizza un'integrazione numerica (metodo di Eulero) per simulare l'evoluzione
 * della tensione del condensatore nel tempo.
 */
class RCOscillatorSimulator {
public:
    /**
     * @brief Costruttore per il simulatore dell'oscillatore RC.
     * @param R_ohms Valore della resistenza in Ohm.
     * @param C_farads Valore della capacità in Farad.
     * @param Vcc_volts Tensione di alimentazione positiva (es. 5.0V).
     * @param Vee_volts Tensione di alimentazione negativa/terra per la scarica (es. 0.0V).
     * @param V_upper_threshold Soglia superiore di tensione per il comparatore.
     * @param V_lower_threshold Soglia inferiore di tensione per il comparatore.
     * @param initial_capacitor_voltage Tensione iniziale del condensatore.
     * @param timestep_seconds Dimensione del passo temporale per la simulazione (es. 0.000001 per 1us).
     */
    RCOscillatorSimulator(double R_ohms, double C_farads, double Vcc_volts, double Vee_volts,
                          double V_upper_threshold, double V_lower_threshold,
                          double initial_capacitor_voltage, double timestep_seconds);

    /**
     * @brief Aggiorna lo stato del circuito per un singolo passo temporale.
     * Questo metodo calcola la nuova tensione del condensatore e lo stato del comparatore.
     */
    void update();

    /**
     * @brief Restituisce la tensione attuale sul condensatore.
     * @return La tensione del condensatore in Volt.
     */
    double getCapacitorVoltage() const { return Vc_; }

    /**
     * @brief Restituisce lo stato di output del comparatore.
     * @return Vcc_volts se il comparatore è alto (carica), Vee_volts se è basso (scarica).
     */
    double getOutputState() const { return comparator_output_state_; }

    // Metodi setter (opzionali per modifiche dinamiche)
    void setResistance(double R_ohms);
    void setCapacitance(double C_farads);
    void setUpperThreshold(double V_upper_threshold);
    void setLowerThreshold(double V_lower_threshold);
    void setTimeStep(double timestep_seconds);

private:
    double R_;                      ///< Resistenza in Ohm.
    double C_;                      ///< Capacità in Farad.
    double Vcc_;                    ///< Tensione di alimentazione positiva.
    double Vee_;                    ///< Tensione di alimentazione negativa/terra.
    double V_upper_;                ///< Soglia superiore del comparatore.
    double V_lower_;                ///< Soglia inferiore del comparatore.
    double Vc_;                     ///< Tensione corrente sul condensatore.
    double dt_;                     ///< Passo temporale della simulazione.

    bool charging_direction_;       ///< Vero se il condensatore si sta caricando verso Vcc, falso se scaricando verso Vee.
    double comparator_output_state_; ///< Stato attuale dell'output del comparatore (Vcc o Vee).
};

#endif // RC_OSCILLATOR_SIMULATOR_H
