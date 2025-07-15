// circuit_solver/MnaSolver.h
#ifndef MNA_SOLVER_H
#define MNA_SOLVER_H

#include <vector>
#include <string>
#include <memory> // Per std::shared_ptr
#include <functional> // Per std::function
#include <iostream> // Per debug

#include "Circuit.h"
#include "components/Component.h"
#include "components/VoltageSource.h"
#include "components/Resistor.h"
#include "components/Capacitor.h"
#include "components/Inductor.h"
#include "components/Diode.h"
#include "components/MOSFET.h"
#include "components/Triode.h"
#include "components/Pentode.h"
#include "components/RectifierTube.h"
#include "components/LED.h"
#include "components/LDR.h"
#include "components/JFET.h"
#include "components/BJT.h"
#include "components/SchottkyDiode.h"
#include "components/ZenerDiode.h"
#include "components/SpeakerDriver.h"
#include "components/ClosedBoxCabinet.h"
#include "components/BassReflexCabinet.h"
#include "components/Splitter.h"

class MnaSolver {
private:
    Circuit& circuit;
    int numNodes;
    int numVoltageSources;
    int numSplitterAuxVars; // Numero di variabili ausiliarie aggiunte dagli splitter
    int numTotalEquations;

    // Liste di componenti classificate per un accesso più rapido
    std::vector<std::shared_ptr<Component>> linearComponents;
    std::vector<std::shared_ptr<Component>> dynamicComponents;
    std::vector<std::shared_ptr<Component>> nonlinearComponents;

    // Metodo per assegnare gli indici delle variabili ausiliarie
    void assignAuxiliaryIndices();
    // Metodo per classificare i componenti
    void classifyComponents();

    // Funzione che rappresenta il sistema di equazioni F(x) = 0 per la risoluzione non lineare
    // x: vettore delle incognite (tensioni ai nodi, correnti Vs, correnti splitter)
    // dt: passo temporale
    // prev_solution: soluzione del passo temporale precedente
    // time: tempo attuale
    // F: vettore dei residui delle equazioni
    void buildSystemEquations(const std::vector<double>& x, double dt, const std::vector<double>& prev_solution, double time,
                              std::vector<double>& F);

    // Implementazione semplificata di Newton-Raphson per fsolve
    // (Per un uso reale, si userebbe una libreria numerica robusta)
    std::vector<double> solveNonlinearSystem(std::vector<double> initial_guess, double dt, const std::vector<double>& prev_solution, double time);
    
    // Calcola la matrice Jacobiana (approssimazione numerica)
    std::vector<std::vector<double>> calculateJacobian(const std::vector<double>& x, double dt, const std::vector<double>& prev_solution, double time, double h = 1e-6);

    // Aggiorna lo stato dei componenti dinamici
    void updateDynamicComponentStates(const std::vector<double>& current_solution, const std::vector<double>& prev_solution, double dt);

public:
    MnaSolver(Circuit& circuit);

    // Esegue una simulazione transitoria del circuito
    std::pair<std::vector<double>, std::vector<std::vector<double>>> simulateTransient(double start_time, double end_time, double time_step);
};

#endif // MNA_SOLVER_H
