// main.cpp
#include <iostream>
#include <vector>
#include <string>
#include <memory> // Per std::shared_ptr
#include <cmath>  // Per M_PI, std::sin, std::exp
#include <iomanip> // Per std::fixed, std::setprecision
#include <algorithm> // Per std::min, std::max

// Includi le classi del solutore e dei componenti
#include "circuit_solver/Circuit.h"
#include "circuit_solver/MnaSolver.h"
#include "circuit_solver/subcircuits/Subcircuit.h" // Per la gestione dei sottocircuiti

// Includi tutti i componenti
#include "components/Component.h" // Base class
#include "components/Resistor.h"
#include "components/Capacitor.h"
#include "components/Inductor.h"
#include "components/VoltageSource.h"
#include "components/CurrentSource.h"
#include "components/Splitter.h"
#include "components/Diode.h"
#include "components/SchottkyDiode.h"
#include "components/ZenerDiode.h"
#include "components/MOSFET.h"
#include "components/Triode.h"
#include "components/Pentode.h"
#include "components/RectifierTube.h"
#include "components/LED.h"
#include "components/JFET.h"
#include "components/BJT.h"
#include "components/LDR.h"
#include "components/SpeakerDriver.h"
#include "components/ClosedBoxCabinet.h"
#include "components/BassReflexCabinet.h"

// --- CLASSI DI SORGENTI PERSONALIZZATE ---
// Queste classi derivano da VoltageSource o CurrentSource per creare segnali specifici.

// Sorgente di tensione sinusoidale
class SineVoltageSource : public VoltageSource {
private:
    double amplitude;
    double frequency;
    double phase; // in radianti

public:
    SineVoltageSource(const std::string& name, const std::string& node_plus, const std::string& node_minus,
                      double amplitude, double frequency, double phase_deg = 0.0)
        : VoltageSource(name, node_plus, node_minus, 0.0, 0.0), // Valore iniziale e tolleranza non usati per segnali variabili
          amplitude(amplitude), frequency(frequency), phase(phase_deg * M_PI / 180.0) {}

    double getVoltage(double time) const override {
        return amplitude * std::sin(2.0 * M_PI * frequency * time + phase);
    }
};

// Sorgente di tensione a onda quadra
class SquareVoltageSource : public VoltageSource {
private:
    double amplitude;
    double frequency;
    double duty_cycle; // Tra 0.0 e 1.0

public:
    SquareVoltageSource(const std::string& name, const std::string& node_plus, const std::string& node_minus,
                        double amplitude, double frequency, double duty_cycle = 0.5)
        : VoltageSource(name, node_plus, node_minus, 0.0, 0.0),
          amplitude(amplitude), frequency(frequency), duty_cycle(duty_cycle) {}

    double getVoltage(double time) const override {
        double period = 1.0 / frequency;
        double t_mod_period = std::fmod(time, period);
        if (t_mod_period < period * duty_cycle) {
            return amplitude;
        } else {
            return 0.0;
        }
    }
};

// Sorgente di corrente impulsiva (impulso singolo)
class PulseCurrentSource : public CurrentSource {
private:
    double amplitude;
    double start_time;
    double duration;

public:
    PulseCurrentSource(const std::string& name, const std::string& node1, const std::string& node2,
                       double amplitude, double start_t, double dur)
        : CurrentSource(name, node1, node2, 0.0, 0.0),
          amplitude(amplitude), start_time(start_t), duration(dur) {}

    double getCurrent(double time) const override {
        if (time >= start_time && time < (start_time + duration)) {
            return amplitude;
        } else {
            return 0.0;
        }
    }
};

// --- FUNZIONI PER LA CREAZIONE DI SOTTOCIRCUITI ---
// Queste funzioni definiscono dei sottocircuiti riutilizzabili.

// Crea una definizione di filtro passa-basso RC
Subcircuit createRcLowPassFilter(const std::string& name, double R_val, double C_val) {
    // Definisce un sottocircuito con tre porte: "input", "output" e "0" (ground)
    Subcircuit rc_filter_def(name, "input", "output", "0");

    // Aggiungi i componenti interni al circuito interno del sottocircuito
    // Nota: "filter_node" è un nodo INTERNO al sottocircuito.
    rc_filter_def.addInternalComponent(std::make_shared<Resistor>(name + "_R", "input", "filter_node", R_val));
    rc_filter_def.addInternalComponent(std::make_shared<Capacitor>(name + "_C", "filter_node", "0", C_val));
    
    // Per mappare il nodo interno "filter_node" alla porta esterna "output",
    // ci affidiamo alla logica di addSubcircuitInstance che userà la mappa di connessione.
    // L'importante è che "output" sia tra i portNames e che "filter_node" sia il nodo interno desiderato.
    // In questo caso, il nome della porta "output" si riferirà al nodo interno "filter_node".
    return rc_filter_def;
}

// Crea una definizione di stadio amplificatore a Triodo (Common Cathode)
Subcircuit createTriodeStage(const std::string& name, double R_plate, double R_cathode, double C_cathode,
                             double triode_K, double triode_Mu, double triode_Ex, double triode_Vg_off) {
    Subcircuit triode_stage_def(name, "input", "output", "V_B+", "0"); // Porte: Input, Output, Alimentazione, Ground

    // Componenti interni
    // Resistore di carico di placca (Plate Resistor)
    triode_stage_def.addInternalComponent(std::make_shared<Resistor>(name + "_Rp", "V_B+", "plate_node", R_plate));
    // Triodo
    triode_stage_def.addInternalComponent(std::make_shared<Triode>(name + "_V1", "plate_node", "input", "cathode_node",
                                                                  triode_K, triode_Mu, triode_Ex, triode_Vg_off));
    // Resistore di catodo (Cathode Resistor)
    triode_stage_def.addInternalComponent(std::make_shared<Resistor>(name + "_Rk", "cathode_node", "0", R_cathode));
    // Condensatore di bypass del catodo (Cathode Bypass Capacitor)
    triode_stage_def.addInternalComponent(std::make_shared<Capacitor>(name + "_Ck", "cathode_node", "0", C_cathode));
    
    // Connessione interna: il nodo "plate_node" è l'output
    // Questo è gestito dalla mappatura delle porte nell'istanza.
    // La porta "output" del sottocircuito si riferirà a "plate_node" internamente.
    return triode_stage_def;
}


int main() {
    std::cout << "--- Simulatore di Circuiti Analogici MNA in C++ ---" << std::endl;
    std::cout << "Esempi di circuiti complessi con vari componenti e sottocircuiti." << std::endl;

    Circuit main_circuit;

    // --- SEZIONE: DEFINIZIONE E UTILIZZO DI SOTTOCIRCUITI ---

    // 1. Definizione di un Filtro Passa-Basso RC riutilizzabile
    std::cout << "\n--- Definizione Sottocircuito: Filtro Passa-Basso RC ---" << std::endl;
    Subcircuit rc_filter_template = createRcLowPassFilter("RC_LPF", 1000.0, 100e-9); // R=1k, C=100nF

    // 2. Istanza 1 del Filtro RC
    std::cout << "\n--- Istanza 1: RC_LPF applicato a un segnale sinusoidale ---" << std::endl;
    ma
