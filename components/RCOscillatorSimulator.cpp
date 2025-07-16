// File: RCOscillatorSimulator.cpp
#include "RCOscillatorSimulator.h"
#include <algorithm> // Per std::max e std::min

/**
 * @brief Costruttore per il simulatore dell'oscillatore RC.
 * @param R_ohms Valore della resistenza in Ohm.
 * @param C_farads Valore della capacità in Farad.
 * @param Vcc_volts Tensione di alimentazione positiva.
 * @param Vee_volts Tensione di alimentazione negativa/terra.
 * @param V_upper_threshold Soglia superiore di tensione per il comparatore.
 * @param V_lower_threshold Soglia inferiore di tensione per il comparatore.
 * @param initial_capacitor_voltage Tensione iniziale del condensatore.
 * @param timestep_seconds Dimensione del passo temporale per la simulazione.
 */
RCOscillatorSimulator::RCOscillatorSimulator(double R_ohms, double C_farads, double Vcc_volts, double Vee_volts,
                                             double V_upper_threshold, double V_lower_threshold,
                                             double initial_capacitor_voltage, double timestep_seconds)
    : R_(R_ohms),
      C_(C_farads),
      Vcc_(Vcc_volts),
      Vee_(Vee_volts),
      V_upper_(V_upper_threshold),
      V_lower_(V_lower_threshold),
      Vc_(initial_capacitor_voltage),
      dt_(timestep_seconds)
{
    // Validazione e clamping dei parametri
    if (R_ <= 0) R_ = 1.0; // Evita divisione per zero
    if (C_ <= 0) C_ = 1e-9; // Capacità minima
    if (dt_ <= 0) dt_ = 1e-6; // Passo temporale minimo

    // Assicurati che le soglie siano nell'intervallo [Vee, Vcc]
    V_lower_ = std::max(Vee_, std::min(Vcc_, V_lower_));
    V_upper_ = std::max(Vee_, std::min(Vcc_, V_upper_));
    if (V_lower_ >= V_upper_) { // Assicurati che la soglia inferiore sia minore della superiore
        V_upper_ = V_lower_ + 0.1; // Aggiungi un piccolo offset se sono uguali o invertite
    }

    // Inizializza lo stato di carica e l'output del comparatore
    // Se la tensione iniziale è sotto la soglia inferiore, inizia a caricare.
    // Altrimenti, inizia a scaricare.
    if (Vc_ <= V_lower_) {
        charging_direction_ = true; // Inizia a caricare
        comparator_output_state_ = Vcc_; // Comparatore alto
    } else {
        charging_direction_ = false; // Inizia a scaricare
        comparator_output_state_ = Vee_; // Comparatore basso
    }
}

/**
 * @brief Aggiorna lo stato del circuito per un singolo passo temporale.
 * Questo metodo calcola la nuova tensione del condensatore e lo stato del comparatore.
 */
void RCOscillatorSimulator::update() {
    // Determina la tensione di riferimento verso cui il condensatore si sta caricando/scaricando
    double target_voltage = charging_direction_ ? Vcc_ : Vee_;

    // Calcola la corrente attraverso il resistore (legge di Ohm)
    // I = (V_sorgente - V_condensatore) / R
    double current_through_resistor = (target_voltage - Vc_) / R_;

    // Calcola la variazione di tensione sul condensatore (legge del condensatore)
    // dV/dt = I/C  =>  dV = (I/C) * dt
    double delta_Vc = (current_through_resistor / C_) * dt_;

    // Aggiorna la tensione del condensatore
    Vc_ += delta_Vc;

    // Gestione delle soglie e cambio di stato
    if (charging_direction_) {
        // Se si sta caricando e si raggiunge o si supera la soglia superiore
        if (Vc_ >= V_upper_) {
            charging_direction_ = false; // Inverti per iniziare a scaricare
            comparator_output_state_ = Vee_; // L'output del comparatore va basso
            Vc_ = V_upper_; // Clampa la tensione per evitare overshoot numerico
        }
    } else { // charging_direction_ è falso, quindi si sta scaricando
        // Se si sta scaricando e si raggiunge o si scende sotto la soglia inferiore
        if (Vc_ <= V_lower_) {
            charging_direction_ = true; // Inverti per iniziare a caricare
            comparator_output_state_ = Vcc_; // L'output del comparatore va alto
            Vc_ = V_lower_; // Clampa la tensione per evitare undershoot numerico
        }
    }

    // Assicurati che la tensione del condensatore rimanga all'interno dei limiti di alimentazione
    Vc_ = std::max(Vee_, std::min(Vcc_, Vc_));
}

// Implementazioni dei metodi setter (opzionali per modifiche dinamiche)
void RCOscillatorSimulator::setResistance(double R_ohms) {
    if (R_ohms > 0) R_ = R_ohms;
}

void RCOscillatorSimulator::setCapacitance(double C_farads) {
    if (C_farads > 0) C_ = C_farads;
}

void RCOscillatorSimulator::setUpperThreshold(double V_upper_threshold) {
    V_upper_ = V_upper_threshold;
    if (V_upper_ <= V_lower_) V_upper_ = V_lower_ + 0.1; // Maintain hysteresis
}

void RCOscillatorSimulator::setLowerThreshold(double V_lower_threshold) {
    V_lower_ = V_lower_threshold;
    if (V_lower_ >= V_upper_) V_lower_ = V_upper_ - 0.1; // Maintain hysteresis
}

void RCOscillatorSimulator::setTimeStep(double timestep_seconds) {
    if (timestep_seconds > 0) dt_ = timestep_seconds;
}
