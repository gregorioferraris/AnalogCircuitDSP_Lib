// components/MOSFET.cpp
#include "MOSFET.h"
#include <iostream> // Per i messaggi di errore e output di debug
#include <algorithm> // Per std::max, std::min
#include <cmath>     // Per std::sqrt, std::abs

// Definisci un piccolo epsilon per la stabilità numerica, specialmente per le derivate vicino ai bordi
const double EPSILON = 1e-9; // Piccola differenza di tensione per la stabilità numerica

// Costruttore per il componente MOSFET.
MOSFET::MOSFET(const std::string& name,
               const std::vector<std::string>& node_names_str, // Ora accetta vector di stringhe
               Type type,
               double Kp, double Vt)
    : Component(name, node_names_str), // Chiama il costruttore della classe base con i nomi dei nodi
      type(type), Kp_(Kp), Vt_(Vt)
{
    // I node_ids_ saranno assegnati dal Circuit.
    // Assumiamo che node_names_str contenga [drain_node, gate_node, source_node, bulk_node]
    // nell'ordine specificato nel costruttore di Component.

    // Validazione di base per i parametri
    if (Kp_ <= 0) {
        std::cerr << "Warning: MOSFET " << name << " ha un parametro Kp non positivo. Impostato a 1e-5 A/V^2." << std::endl;
        Kp_ = 1e-5; // Valore di default
    }

    std::cout << "MOSFET " << name << " inizializzato (Tipo: " << (type == N_CHANNEL ? "NMOS" : "PMOS")
              << ", Kp=" << Kp_ << ", Vt=" << Vt_ << "V)." << std::endl;
}

// Calcola la tensione di soglia (Vth) senza considerare l'effetto corpo (per semplicità Level 1)
// Il tuo file originale MOSFET.cpp aveva l'effetto corpo. Per allineare con la semplificazione
// del tuo MnaSolver.cpp che non usa Vbs, lo rimuovo qui per ora.
// Se vuoi l'effetto corpo, dovrai assicurarti che Vbs sia passato correttamente e che
// MnaSolver.cpp lo utilizzi per calcolare Vth.
double MOSFET::calculateThresholdVoltage(double Vsb) const {
    // Per SPICE Level 1 semplificato senza effetto corpo, Vth è semplicemente Vt0
    // Se vuoi l'effetto corpo, ripristina la logica dal tuo file originale.
    // double phi = 0.6; // Esempio di phi, dovrebbe essere un parametro del modello
    // double gamma = 0.0; // Esempio di gamma, dovrebbe essere un parametro del modello
    // if (type == N_CHANNEL) {
    //     double Vsb_clamped = std::max(0.0, Vsb);
    //     return Vt_ + gamma * (std::sqrt(2 * phi + Vsb_clamped) - std::sqrt(2 * phi));
    // } else { // P_CHANNEL
    //     double Vbs_abs = std::abs(Vsb);
    //     double Vbs_clamped = std::max(0.0, Vbs_abs);
    //     return Vt_ - gamma * (std::sqrt(2 * phi + Vbs_clamped) - std::sqrt(2 * phi));
    // }
    return Vt_;
}

// Calcola la corrente di drain (Id) basata sulla regione operativa (SPICE Level 1).
double MOSFET::calculateDrainCurrent(double Vgs, double Vds, double Vth) const {
    double id = 0.0;
    double Vov = Vgs - Vth; // Tensione di overdrive

    if (type == N_CHANNEL) {
        if (Vgs <= Vth + EPSILON) { // Regione di cutoff
            id = 0.0;
        } else if (Vds < Vov - EPSILON) { // Regione di triodo (lineare)
            id = Kp_ * (Vov * Vds - 0.5 * Vds * Vds); // Senza lambda per semplicità Level 1
        } else { // Regione di saturazione (Vds >= Vov)
            id = 0.5 * Kp_ * Vov * Vov; // Senza lambda per semplicità Level 1
        }
    } else { // P_CHANNEL
        // Per PMOS, le tensioni sono tipicamente negative.
        // Convertiamo in valori assoluti per usare le stesse equazioni, poi invertiamo la corrente.
        double abs_Vgs = std::abs(Vgs);
        double abs_Vds = std::abs(Vds);
        double abs_Vth = std::abs(Vth); // Vth è negativo per PMOS

        if (Vgs >= Vth - EPSILON) { // Regione di cutoff (Vgs meno negativo di Vth)
            id = 0.0;
        } else if (Vds > Vov + EPSILON) { // Regione di triodo (Vds meno negativo di Vds_sat)
            id = Kp_ * (Vov * Vds - 0.5 * Vds * Vds);
        } else { // Regione di saturazione (Vds <= Vov)
            id = 0.5 * Kp_ * Vov * Vov;
        }
        id = -id; // La corrente fluisce da Source a Drain per PMOS, quindi invertiamo il segno
    }
    return id;
}

// Calcola la conduttanza di uscita (gds = d(Id)/d(Vds)).
double MOSFET::calculateGDS(double Vgs, double Vds, double Vth) const {
    double gds = 0.0;
    double Vov = Vgs - Vth;

    if (type == N_CHANNEL) {
        if (Vgs <= Vth + EPSILON) { // Cutoff
            gds = 0.0;
        } else if (Vds < Vov - EPSILON) { // Triodo
            gds = Kp_ * (Vov - Vds);
        } else { // Saturation
            gds = 0.0; // Per Level 1 senza lambda
        }
    } else { // P_CHANNEL
        if (Vgs >= Vth - EPSILON) { // Cutoff
            gds = 0.0;
        } else if (Vds > Vov + EPSILON) { // Triodo
            gds = -Kp_ * (Vov - Vds); // Invertito per PMOS
        } else { // Saturation
            gds = 0.0; // Per Level 1 senza lambda
        }
    }
    return gds;
}

// Calcola la transconduttanza (gm = d(Id)/d(Vgs)).
double MOSFET::calculateGM(double Vgs, double Vds, double Vth) const {
    double gm = 0.0;
    double Vov = Vgs - Vth;

    if (type == N_CHANNEL) {
        if (Vgs <= Vth + EPSILON) { // Cutoff
            gm = 0.0;
        } else if (Vds < Vov - EPSILON) { // Triodo
            gm = Kp_ * Vds;
        } else { // Saturation
            gm = Kp_ * Vov;
        }
    } else { // PMOS
        if (Vgs >= Vth - EPSILON) { // Cutoff
            gm = 0.0;
        } else if (Vds > Vov + EPSILON) { // Triodo
            gm = -Kp_ * Vds; // Invertito per PMOS
        } else { // Saturation
            gm = -Kp_ * Vov; // Invertito per PMOS
        }
    }
    return gm;
}

// Calcola la transconduttanza di bulk (gmb = d(Id)/d(Vbs)).
// Per semplicità, in questo modello Level 1 senza effetto corpo, gmb è 0.
double MOSFET::calculateGMB(double Vgs, double Vds, double Vbs, double Vth) const {
    return 0.0; // Implementazione semplificata senza effetto corpo
}

// Applica gli "stamps" del modello MOSFET alla matrice MNA (A) e al vettore (B).
void MOSFET::getStamps(
    int num_total_equations, double dt,
    const std::vector<double>& x,
    const std::vector<double>& prev_solution,
    double time,
    std::vector<std::vector<double>>& A,
    std::vector<double>& B
) {
    // Ottieni gli indici globali per i nodi del MOSFET.
    // Assumiamo che node_ids_ contenga [idx_D, idx_G, idx_S, idx_B]
    // nell'ordine: Drain (0), Gate (1), Source (2), Bulk (3)
    int idx_D = node_ids_[0];
    int idx_G = node_ids_[1];
    int idx_S = node_ids_[2];
    int idx_B = node_ids_[3];

    // Ottieni le tensioni correnti dei nodi dal vettore `x`
    // Il nodo 0 è il ground, quindi la sua tensione è 0.0
    double V_D_curr = (idx_D == 0) ? 0.0 : x[idx_D];
    double V_G_curr = (idx_G == 0) ? 0.0 : x[idx_G];
    double V_S_curr = (idx_S == 0) ? 0.0 : x[idx_S];
    double V_B_curr = (idx_B == 0) ? 0.0 : x[idx_B];

    // Calcola le tensioni terminali
    double Vgs_curr = V_G_curr - V_S_curr;
    double Vds_curr = V_D_curr - V_S_curr;
    double Vbs_curr = V_B_curr - V_S_curr;
    double Vsb_curr = V_S_curr - V_B_curr; // Source-Bulk voltage per il calcolo di Vth (se usato)

    // Calcola la tensione di soglia (Vth)
    double Vth_curr = calculateThresholdVoltage(Vsb_curr);

    // Calcola la corrente di drain e le conduttanze a piccolo segnale
    double Id_curr = calculateDrainCurrent(Vgs_curr, Vds_curr, Vth_curr);
    double gds_val = calculateGDS(Vgs_curr, Vds_curr, Vth_curr);
    double gm_val = calculateGM(Vgs_curr, Vds_curr, Vth_curr);
    double gmb_val = calculateGMB(Vgs_curr, Vds_curr, Vbs_curr, Vth_curr); // Sarà 0.0 nel modello semplificato

    // Applica gli stamps alla matrice MNA A e al vettore B
    // Il modello del compagno per una sorgente di corrente non lineare I_D(V_GS, V_DS, V_BS) è:
    // I_D = I_D(V_0) + gm * (V_GS - V_GS0) + gds * (V_DS - V_DS0) + gmb * (V_BS - V_BS0)
    // Riorganizzando per MNA:
    // I_D - gm * V_GS - gds * V_DS - gmb * V_BS = I_D(V_0) - gm * V_GS0 - gds * V_DS0 - gmb * V_BS0
    // Il lato destro è la sorgente di corrente equivalente I_eq.

    // Contributi alla matrice A (conduttanze)
    // KCL al nodo Drain (idx_D)
    A[idx_D][idx_D] += gds_val;
    A[idx_D][idx_G] += gm_val;
    A[idx_D][idx_B] += gmb_val;
    A[idx_D][idx_S] -= (gds_val + gm_val + gmb_val);

    // KCL al nodo Source (idx_S)
    A[idx_S][idx_D] -= gds_val;
    A[idx_S][idx_G] -= gm_val;
    A[idx_S][idx_B] -= gmb_val;
    A[idx_S][idx_S] += (gds_val + gm_val + gmb_val);

    // Calcola la sorgente di corrente equivalente per il vettore RHS B
    double I_eq = Id_curr - gds_val * Vds_curr - gm_val * Vgs_curr - gmb_val * Vbs_curr;

    // Aggiungi I_eq al RHS dell'equazione KCL del Drain
    // La corrente Id_curr è la corrente che *esce* dal Drain per NMOS, o *entra* per PMOS.
    // Quindi, se Id_curr è la corrente che fluisce da D a S (NMOS), la sottraiamo dal nodo D.
    // Se Id_curr è la corrente che fluisce da S a D (PMOS), la aggiungiamo al nodo D.
    // Poiché calculateDrainCurrent() restituisce un valore positivo per NMOS (D->S)
    // e un valore negativo per PMOS (S->D), possiamo semplicemente sottrarre I_eq dal Drain
    // e aggiungerlo al Source.
    B[idx_D] -= I_eq; // Corrente esce dal Drain
    B[idx_S] += I_eq; // Corrente entra nel Source

    // Nota: I nodi Gate e Bulk sono tipicamente ad alta impedenza (nessun contributo di corrente)
    // a meno che non siano modellati la corrente di gate o i diodi di giunzione,
    // che non fanno parte del modello Level 1 di base.
}

// updateState non è necessario per il MOSFET nel modello Level 1,
// poiché non ha variabili di stato interne che cambiano nel tempo
// in base alla soluzione precedente (come condensatori/induttori).
// Il suo comportamento è puramente non lineare e dipende dalle tensioni attuali.
void MOSFET::updateState(const std::vector<double>& current_solution,
                         const std::vector<double>& prev_solution,
                         double dt) {
    // Non fa nulla per il MOSFET Level 1, in quanto non ha stato interno dinamico.
}
