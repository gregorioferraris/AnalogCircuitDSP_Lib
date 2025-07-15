// circuit_solver/MnaSolver.cpp
#include "MnaSolver.h"
#include <iomanip> // Per std::fixed, std::setprecision
#include <cmath>   // Per std::abs, std::sqrt, std::nan

// Funzione helper per l'inversione di matrice (molto basilare, non per grandi matrici)
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& A) {
    int n = static_cast<int>(A.size());
    if (n == 0 || A[0].size() != n) {
        throw std::invalid_argument("Matrix must be square and non-empty.");
    }

    std::vector<std::vector<double>> identity(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> augmented_A = A;

    for (int i = 0; i < n; ++i) {
        identity[i][i] = 1.0;
        augmented_A[i].insert(augmented_A[i].end(), identity[i].begin(), identity[i].end());
    }

    // Eliminazione di Gauss-Jordan
    for (int i = 0; i < n; ++i) {
        // Trova il pivot
        int pivot_row = i;
        for (int j = i + 1; j < n; ++j) {
            if (std::abs(augmented_A[j][i]) > std::abs(augmented_A[pivot_row][i])) {
                pivot_row = j;
            }
        }
        std::swap(augmented_A[i], augmented_A[pivot_row]);

        double pivot = augmented_A[i][i];
        if (std::abs(pivot) < 1e-12) { // Matrice singolare o quasi singolare
            throw std::runtime_error("Matrix is singular or ill-conditioned.");
        }

        for (int j = i; j < 2 * n; ++j) {
            augmented_A[i][j] /= pivot;
        }

        for (int row = 0; row < n; ++row) {
            if (row != i) {
                double factor = augmented_A[row][i];
                for (int col = i; col < 2 * n; ++col) {
                    augmented_A[row][col] -= factor * augmented_A[i][col];
                }
            }
        }
    }

    std::vector<std::vector<double>> inv_A(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inv_A[i][j] = augmented_A[i][j + n];
        }
    }
    return inv_A;
}

// Funzione helper per la moltiplicazione matrice-vettore
std::vector<double> multiply(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    int rows = static_cast<int>(A.size());
    if (rows == 0) return {};
    int cols = static_cast<int>(A[0].size());
    if (cols != x.size()) {
        throw std::invalid_argument("Matrix columns must match vector rows.");
    }

    std::vector<double> result(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i] += A[i][j] * x[j];
        }
    }
    return result;
}

MnaSolver::MnaSolver(Circuit& circuit)
    : circuit(circuit),
      numNodes(circuit.getNumNodes()),
      numVoltageSources(static_cast<int>(circuit.getVoltageSources().size())),
      numSplitterAuxVars(0) // Inizializzato a 0, calcolato in assignAuxiliaryIndices
{
    // Calcola il numero di variabili ausiliarie per gli splitter
    for (const auto& splitter : circuit.getSplitters()) {
        numSplitterAuxVars += splitter->num_outputs;
    }

    numTotalEquations = numNodes + numVoltageSources + numSplitterAuxVars;

    assignAuxiliaryIndices();
    classifyComponents();

    std::cout << "MnaSolver inizializzato. Equazioni totali: " << numTotalEquations << std::endl;
}

void MnaSolver::assignAuxiliaryIndices() {
    int currentAuxIdx = numNodes; // Inizia dopo gli ID dei nodi

    // Assegna indici alle correnti delle sorgenti di tensione
    for (const auto& vs : circuit.getVoltageSources()) {
        vs->set_current_index(currentAuxIdx);
        currentAuxIdx++;
    }

    // Assegna indici alle correnti ausiliarie degli splitter
    for (const auto& splitter : circuit.getSplitters()) {
        std::vector<int> indices;
        for (int i = 0; i < splitter->num_outputs; ++i) {
            indices.push_back(currentAuxIdx);
            currentAuxIdx++;
        }
        splitter->_set_output_current_indices(indices);
    }
}

void MnaSolver::classifyComponents() {
    for (const auto& comp : circuit.getComponents()) {
        if (std::dynamic_pointer_cast<Diode>(comp) ||
            std::dynamic_pointer_cast<SchottkyDiode>(comp) ||
            std::dynamic_pointer_cast<ZenerDiode>(comp) ||
            std::dynamic_pointer_cast<MOSFET>(comp) ||
            std::dynamic_pointer_cast<Triode>(comp) ||
            std::dynamic_pointer_cast<Pentode>(comp) ||
            std::dynamic_pointer_cast<RectifierTube>(comp) ||
            std::dynamic_pointer_cast<LED>(comp) ||
            std::dynamic_pointer_cast<JFET>(comp) ||
            std::dynamic_pointer_cast<BJT>(comp)) {
            nonlinearComponents.push_back(comp);
        } else if (std::dynamic_pointer_cast<Capacitor>(comp) ||
                   std::dynamic_pointer_cast<Inductor>(comp) ||
                   std::dynamic_pointer_cast<SpeakerDriver>(comp) ||
                   std::dynamic_pointer_cast<ClosedBoxCabinet>(comp) ||
                   std::dynamic_pointer_cast<BassReflexCabinet>(comp) ||
                   std::dynamic_pointer_cast<LDR>(comp)) {
            dynamicComponents.push_back(comp);
        } else { // Resistor, VoltageSource, CurrentSource, Splitter
            linearComponents.push_back(comp);
        }
    }
}

void MnaSolver::buildSystemEquations(const std::vector<double>& x, double dt, const std::vector<double>& prev_solution, double time,
                                     std::vector<double>& F) {
    // Inizializza la matrice MNA (A) e il vettore RHS (B)
    std::vector<std::vector<double>> A(numTotalEquations, std::vector<double>(numTotalEquations, 0.0));
    std::vector<double> B(numTotalEquations, 0.0);

    // --- Contributi dei componenti lineari e dinamici ---
    for (const auto& comp : linearComponents) {
        comp->getStamps(numTotalEquations, dt, x, prev_solution, time, A, B);
    }
    for (const auto& comp : dynamicComponents) {
        comp->getStamps(numTotalEquations, dt, x, prev_solution, time, A, B);
    }

    // --- Contributi dei componenti non lineari (gestiti implicitamente) ---
    std::vector<double> I_nonlinear(numTotalEquations, 0.0);

    for (const auto& comp : nonlinearComponents) {
        const std::vector<int>& node_ids = comp->getNodeIds();
        
        if (std::shared_ptr<Diode> diode = std::dynamic_pointer_cast<Diode>(comp)) {
            double Vd = x[node_ids[0]] - x[node_ids[1]];
            double current = diode->calculateCurrent(Vd);
            if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += current;
            if (node_ids[1] != 0) I_nonlinear[node_ids[1]] -= current;
        } else if (std::shared_ptr<SchottkyDiode> schottky = std::dynamic_pointer_cast<SchottkyDiode>(comp)) {
            double Vd = x[node_ids[0]] - x[node_ids[1]];
            double current = schottky->calculateCurrent(Vd);
            if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += current;
            if (node_ids[1] != 0) I_nonlinear[node_ids[1]] -= current;
        } else if (std::shared_ptr<ZenerDiode> zener = std::dynamic_pointer_cast<ZenerDiode>(comp)) {
            double Vd = x[node_ids[0]] - x[node_ids[1]];
            double current = zener->calculateCurrent(Vd);
            if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += current;
            if (node_ids[1] != 0) I_nonlinear[node_ids[1]] -= current;
        } else if (std::shared_ptr<MOSFET> mosfet = std::dynamic_pointer_cast<MOSFET>(comp)) {
            double Vds = x[node_ids[0]] - x[node_ids[2]]; // Drain - Source
            double Vgs = x[node_ids[1]] - x[node_ids[2]]; // Gate - Source
            double Id = mosfet->calculateDrainCurrent(Vgs, Vds);
            if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += Id; // Corrente entra nel Drain
            if (node_ids[2] != 0) I_nonlinear[node_ids[2]] -= Id; // Corrente esce dal Source
        } else if (std::shared_ptr<Triode> triode = std::dynamic_pointer_cast<Triode>(comp)) {
            double Vpk = x[node_ids[0]] - x[node_ids[2]]; // Plate - Cathode
            double Vgk = x[node_ids[1]] - x[node_ids[2]]; // Grid - Cathode
            double Ip = triode->calculatePlateCurrent(Vgk, Vpk);
            if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += Ip; // Corrente entra nella Placca
            if (node_ids[2] != 0) I_nonlinear[node_ids[2]] -= Ip; // Corrente esce dal Catodo
        } else if (std::shared_ptr<Pentode> pentode = std::dynamic_pointer_cast<Pentode>(comp)) {
            // Assumi che node_ids sia nell'ordine (plate, grid, screen_grid, suppressor_grid, cathode)
            double Vpk = x[node_ids[pentode->pin_map["plate"]]] - x[node_ids[pentode->pin_map["cathode"]]];
            double Vgk = x[node_ids[pentode->pin_map["grid"]]] - x[node_ids[pentode->pin_map["cathode"]]];
            double Vg2k = x[node_ids[pentode->pin_map["screen_grid"]]] - x[node_ids[pentode->pin_map["cathode"]]];
            double Vg3k = x[node_ids[pentode->pin_map["suppressor_grid"]]] - x[node_ids[pentode->pin_map["cathode"]]];
            double Ip = pentode->calculatePlateCurrent(Vgk, Vpk, Vg2k, Vg3k);
            if (node_ids[pentode->pin_map["plate"]] != 0) I_nonlinear[node_ids[pentode->pin_map["plate"]]] += Ip;
            if (node_ids[pentode->pin_map["cathode"]] != 0) I_nonlinear[node_ids[pentode->pin_map["cathode"]]] -= Ip;
        } else if (std::shared_ptr<RectifierTube> rectifier = std::dynamic_pointer_cast<RectifierTube>(comp)) {
            double Vd = x[node_ids[0]] - x[node_ids[1]]; // Anode - Cathode
            double current = rectifier->calculateCurrent(Vd);
            if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += current;
            if (node_ids[1] != 0) I_nonlinear[node_ids[1]] -= current;
        } else if (std::shared_ptr<LED> led = std::dynamic_pointer_cast<LED>(comp)) {
            double Vd = x[node_ids[0]] - x[node_ids[1]]; // Anode - Cathode
            double current = led->calculateCurrent(Vd);
            if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += current;
            if (node_ids[1] != 0) I_nonlinear[node_ids[1]] -= current;
        } else if (std::shared_ptr<JFET> jfet = std::dynamic_pointer_cast<JFET>(comp)) {
            double Vds = x[node_ids[0]] - x[node_ids[2]]; // Drain - Source
            double Vgs = x[node_ids[1]] - x[node_ids[2]]; // Gate - Source
            double Id = jfet->calculateDrainCurrent(Vgs, Vds);
            if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += Id;
            if (node_ids[2] != 0) I_nonlinear[node_ids[2]] -= Id;
        } else if (std::shared_ptr<BJT> bjt = std::dynamic_pointer_cast<BJT>(comp)) {
            double Vbe = x[node_ids[1]] - x[node_ids[2]]; // Base - Emitter
            double Vbc = x[node_ids[1]] - x[node_ids[0]]; // Base - Collector
            auto currents = bjt->calculateCurrents(Vbe, Vbc);
            double Ic = currents.first;
            double Ib = currents.second;
            double Ie = Ic + Ib; // Corrente di emettitore
            
            if (bjt->type == BJT::NPN) {
                if (node_ids[0] != 0) I_nonlinear[node_ids[0]] += Ic; // Ic entra nel Collector
                if (node_ids[1] != 0) I_nonlinear[node_ids[1]] += Ib; // Ib entra nella Base
                if (node_ids[2] != 0) I_nonlinear[node_ids[2]] -= Ie; // Ie esce dall'Emitter
            } else if (bjt->type == BJT::PNP) {
                if (node_ids[0] != 0) I_nonlinear[node_ids[0]] -= Ic; // Ic esce dal Collector
                if (node_ids[1] != 0) I_nonlinear[node_ids[1]] -= Ib; // Ib esce dalla Base
                if (node_ids[2] != 0) I_nonlinear[node_ids[2]] += Ie; // Ie entra nell'Emitter
            }
        }
    }

    // Costruisci il vettore delle equazioni F(x) = A*x - B - I_nonlinear = 0
    F.assign(numTotalEquations, 0.0); // Assicurati che F abbia la dimensione corretta e sia inizializzato a zero

    for (int i = 0; i < numTotalEquations; ++i) {
        for (int j = 0; j < numTotalEquations; ++j) {
            F[i] += A[i][j] * x[j];
        }
        F[i] -= B[i];
        F[i] -= I_nonlinear[i];
    }

    // Impone V_ground = 0
    F[0] = x[0];
}

std::vector<double> MnaSolver::calculateJacobian(const std::vector<double>& x, double dt, const std::vector<double>& prev_solution, double time, double h) {
    std::vector<std::vector<double>> J(numTotalEquations, std::vector<double>(numTotalEquations, 0.0));
    std::vector<double> F_base(numTotalEquations);
    buildSystemEquations(x, dt, prev_solution, time, F_base);

    for (int j = 0; j < numTotalEquations; ++j) {
        std::vector<double> x_perturbed = x;
        x_perturbed[j] += h;
        std::vector<double> F_perturbed(numTotalEquations);
        buildSystemEquations(x_perturbed, dt, prev_solution, time, F_perturbed);

        for (int i = 0; i < numTotalEquations; ++i) {
            J[i][j] = (F_perturbed[i] - F_base[i]) / h;
        }
    }
    return J;
}

std::vector<double> MnaSolver::solveNonlinearSystem(std::vector<double> initial_guess, double dt, const std::vector<double>& prev_solution, double time) {
    std::vector<double> x_k = initial_guess;
    const int max_iterations = 100;
    const double tolerance = 1e-9;

    for (int k = 0; k < max_iterations; ++k) {
        std::vector<double> F_k(numTotalEquations);
        buildSystemEquations(x_k, dt, prev_solution, time, F_k);

        double norm_F = 0.0;
        for (double val : F_k) {
            norm_F += val * val;
        }
        norm_F = std::sqrt(norm_F);

        if (norm_F < tolerance) {
            return x_k; // Convergenza raggiunta
        }

        std::vector<std::vector<double>> J_k = calculateJacobian(x_k, dt, prev_solution, time);

        // Risolvi J_k * delta_x = -F_k
        try {
            std::vector<std::vector<double>> J_k_inv = inverse(J_k);
            std::vector<double> neg_F_k = F_k;
            for (double& val : neg_F_k) {
                val = -val;
            }
            std::vector<double> delta_x = multiply(J_k_inv, neg_F_k);

            // Aggiorna x_k+1 = x_k + delta_x
            for (int i = 0; i < numTotalEquations; ++i) {
                x_k[i] += delta_x[i];
            }
        } catch (const std::runtime_error& e) {
            std::cerr << "Errore Jacobiano a iterazione " << k << ": " << e.what() << std::endl;
            return initial_guess; // Ritorna il guess iniziale in caso di errore
        }
    }
    std::cerr << "Avviso: Newton-Raphson non converge dopo " << max_iterations << " iterazioni." << std::endl;
    return x_k; // Ritorna l'ultima approssimazione
}


void MnaSolver::updateDynamicComponentStates(const std::vector<double>& current_solution, const std::vector<double>& prev_solution, double dt) {
    for (const auto& comp : dynamicComponents) {
        if (std::shared_ptr<Capacitor> cap = std::dynamic_pointer_cast<Capacitor>(comp)) {
            int node1_id = cap->getNodeIds()[0];
            int node2_id = cap->getNodeIds()[1];
            double V_curr = current_solution[node1_id] - current_solution[node2_id];
            double V_prev = prev_solution[node1_id] - prev_solution[node2_id];
            double I_curr = (2.0 * cap->capacitance / dt) * (V_curr - V_prev) - cap->i_prev;
            cap->updateState(V_curr, I_curr);
        } else if (std::shared_ptr<Inductor> ind = std::dynamic_pointer_cast<Inductor>(comp)) {
            int node1_id = ind->getNodeIds()[0];
            int node2_id = ind->getNodeIds()[1];
            double V_curr = current_solution[node1_id] - current_solution[node2_id];
            double V_prev = prev_solution[node1_id] - prev_solution[node2_id];
            double I_curr = ind->i_prev + (dt / (2.0 * ind->L)) * (V_curr + V_prev);
            ind->updateState(V_curr, I_curr);
        } else if (std::shared_ptr<SpeakerDriver> spk = std::dynamic_pointer_cast<SpeakerDriver>(comp)) {
            // Per SpeakerDriver, i nodi sono (elec_plus, elec_minus, mech_velocity_out)
            // La corrente attraverso Le (i_Le_curr) e la tensione ai suoi capi (v_Le_curr)
            // La tensione ai capi di Cms (v_Cms_curr) e la corrente attraverso (i_Cms_curr)
            // La tensione ai capi di Mms (v_Mms_curr) e la corrente attraverso (i_Mms_curr)
            // Questi dovrebbero essere estratti da current_solution e passati.
            // La modellazione completa richiede nodi interni espliciti per ciascun elemento meccanico equivalente.
            // Per ora, un placeholder, assumendo che i_Le_curr sia la corrente totale attraverso la bobina.
            // Questo è un punto che richiederebbe un'implementazione più dettagliata del modello.
            
            // Esempio semplificato:
            int elec_plus_id = spk->getNodeIds()[0];
            int elec_minus_id = spk->getNodeIds()[1];
            int mech_velocity_out_id = spk->getNodeIds()[2];

            double v_Le_curr_val = current_solution[elec_plus_id] - current_solution[elec_minus_id];
            // Questo è un placeholder. La corrente effettiva attraverso Le andrebbe calcolata dal solutore.
            double i_Le_curr_val = 0.0; // Placeholder, da calcolare
            
            double v_Cms_curr_val = current_solution[mech_velocity_out_id] - current_solution[0];
            double i_Cms_curr_val = 0.0; // Placeholder
            
            double v_Mms_curr_val = current_solution[mech_velocity_out_id] - current_solution[0];
            double i_Mms_curr_val = 0.0; // Placeholder

            spk->updateState(i_Le_curr_val, v_Le_curr_val, v_Cms_curr_val, i_Cms_curr_val, v_Mms_curr_val, i_Mms_curr_val);
        } else if (std::shared_ptr<ClosedBoxCabinet> cb = std::dynamic_pointer_cast<ClosedBoxCabinet>(comp)) {
            int output_pressure_id = cb->getNodeIds()[1]; // Pressione acustica
            double V_curr = current_solution[output_pressure_id] - current_solution[0];
            double V_prev = prev_solution[output_pressure_id] - prev_solution[0];
            double I_curr = (2.0 * cb->C_box / dt) * (V_curr - V_prev) - cb->i_Cbox_prev;
            cb->updateState(V_curr, I_curr);
        } else if (std::shared_ptr<BassReflexCabinet> br = std::dynamic_pointer_cast<BassReflexCabinet>(comp)) {
            int output_pressure_id = br->getNodeIds()[1]; // Pressione acustica
            double V_Cbox_curr = current_solution[output_pressure_id] - current_solution[0];
            double V_Cbox_prev = prev_solution[output_pressure_id] - prev_solution[0];
            double I_Cbox_curr = (2.0 * br->C_box / dt) * (V_Cbox_curr - V_Cbox_prev) - br->i_Cbox_prev;

            // Per L_port e R_port, avresti bisogno di nodi interni per la corrente e la tensione.
            // Qui, assumiamo che la tensione ai capi del condotto sia la stessa del volume della cassa.
            double V_Lport_curr = V_Cbox_curr;
            double V_Lport_prev = V_Cbox_prev;
            double I_Lport_curr = br->i_Lport_prev + (dt / (2.0 * br->L_port)) * (V_Lport_curr + V_Lport_prev);

            br->updateState(V_Cbox_curr, I_Cbox_curr, V_Lport_curr, I_Lport_curr);
        } else if (std::shared_ptr<LDR>(comp)) {
            // LDR non ha uno stato dinamico interno che cambia con la soluzione MNA,
            // la sua resistenza è impostata esternamente. Nessun update_state qui.
        }
    }
}


std::pair<std::vector<double>, std::vector<std::vector<double>>> MnaSolver::simulateTransient(double start_time, double end_time, double time_step) {
    std::vector<double> times;
    for (double t = start_time; t <= end_time + time_step * 0.5; t += time_step) { // Aggiungi 0.5*time_step per includere end_time
        times.push_back(t);
    }
    
    std::vector<double> prev_solution(numTotalEquations, 0.0); // Inizializza a zero
    std::vector<std::vector<double>> solution_history;

    std::cout << "\nInizio simulazione transitoria da " << std::fixed << std::setprecision(6) << start_time << "s a " << end_time << "s con dt=" << time_step << "s." << std::endl;
    std::cout << "Numero totale di equazioni MNA: " << numTotalEquations << std::endl;

    for (size_t i = 0; i < times.size(); ++i) {
        double t = times[i];
        std::vector<double> initial_guess = prev_solution; // Guess iniziale è la soluzione precedente

        try {
            std::vector<double> current_solution = solveNonlinearSystem(initial_guess, time_step, prev_solution, t);
            
            updateDynamicComponentStates(current_solution, prev_solution, time_step);
            
            prev_solution = current_solution;
            solution_history.push_back(current_solution);
            
            // Debug output
            // if (i % 100 == 0) { // Stampa ogni 100 passi per non intasare la console
            //     std::cout << "Tempo: " << std::fixed << std::setprecision(6) << t << "s";
            //     // Esempio di stampa di una tensione di un nodo specifico (se esiste)
            //     auto it_node_map = circuit.getNodeMap().find("output_final");
            //     if (it_node_map != circuit.getNodeMap().end()) {
            //         std::cout << ", V_out: " << current_solution[it_node_map->second] << "V";
            //     }
            //     std::cout << std::endl;
            // }

        } catch (const std::exception& e) {
            std::cerr << "Errore durante la risoluzione a t=" << std::fixed << std::setprecision(6) << t << "s: " << e.what() << std::endl;
            // Aggiungi NaN o l'ultima soluzione valida in caso di errore
            solution_history.push_back(std::vector<double>(numTotalEquations, std::nan("")));
            break; // Interrompi la simulazione in caso di errore critico
        }
    }
    
    // Assicurati che times e solution_history abbiano la stessa dimensione in caso di interruzione
    if (solution_history.size() < times.size()) {
        times.resize(solution_history.size());
    }

    return {times, solution_history};
}
