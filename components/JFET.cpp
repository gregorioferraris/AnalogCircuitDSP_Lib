// components/JFET.cpp
#include "JFET.h"
#include <iostream> // Per i messaggi di errore
#include <limits>   // Per std::numeric_limits
#include <cmath>    // Per std::abs, std::sqrt

// Costruttore per il componente JFET.
JFET::JFET(const std::string& name,
           const std::vector<std::string>& node_names_str, // Ora accetta vector di stringhe
           Type type, double Idss, double Vp)
    : Component(name, node_names_str), // Chiama il costruttore della classe base con i nomi dei nodi
      type(type), Idss_(Idss), Vp_(Vp),
      prev_VGS_(0.0) // Inizializza la variabile di stato interna
{
    // Il costruttore base Component ora gestisce node_names_str.
    // I node_ids_ saranno assegnati dal Circuit.

    // Validazione dei parametri di input
    if (Idss_ <= 0.0) {
        std::cerr << "Warning: JFET " << name << " ha un Idss non positivo. Impostato a 1e-3A." << std::endl;
        this->Idss_ = 1e-3; // Valore tipico di default
    }
    // Per JFET a canale N, Vp dovrebbe essere negativo. Per JFET a canale P, Vp dovrebbe essere positivo.
    if ((type == N_CHANNEL && Vp_ >= 0.0) || (type == P_CHANNEL && Vp_ <= 0.0)) {
        std::cerr << "Warning: JFET " << name << " ha un Vp (tensione di pinch-off) non coerente con il tipo di canale. "
                  << "Impostato a " << (type == N_CHANNEL ? -2.0 : 2.0) << "V." << std::endl;
        this->Vp_ = (type == N_CHANNEL ? -2.0 : 2.0); // Valore tipico di default
    }

    std::cout << "JFET " << name << " inizializzato (Tipo: " << (type == N_CHANNEL ? "N-Canale" : "P-Canale")
              << ", Idss=" << Idss_ << "A, Vp=" << Vp_ << "V)." << std::endl;
}

// Funzione helper per calcolare la corrente di drain (ID) basata su VGS.
double JFET::calculateDrainCurrent(double VGS) const {
    double id = 0.0;
    if (type == N_CHANNEL) {
        if (VGS >= Vp_) { // Regione di cutoff per N-canale
            id = 0.0;
        } else { // Regione di saturazione per N-canale
            double term = (1.0 - VGS / Vp_);
            id = Idss_ * term * term;
        }
    } else { // P_CHANNEL
        if (VGS <= Vp_) { // Regione di cutoff per P-canale (VGS più positivo di Vp)
            id = 0.0;
        } else { // Regione di saturazione per P-canale
            double term = (1.0 - VGS / Vp_);
            id = -Idss_ * term * term; // Corrente in direzione opposta
        }
    }
    return id;
}

// Funzione helper per calcolare la transconduttanza (gm) basata su VGS.
double JFET::calculateTransconductance(double VGS) const {
    double gm = 0.0;
    if (type == N_CHANNEL) {
        if (VGS >= Vp_) { // Regione di cutoff
            gm = 0.0;
        } else { // Regione di saturazione
            gm = -2.0 * Idss_ / Vp_ * (1.0 - VGS / Vp_);
        }
    } else { // P_CHANNEL
        if (VGS <= Vp_) { // Regione di cutoff
            gm = 0.0;
        } else { // Regione di saturazione
            gm = -2.0 * (-Idss_) / Vp_ * (1.0 - VGS / Vp_); // gm per PMOS
        }
    }
    return gm;
}

// Applica gli "stamps" del modello JFET alla matrice MNA (A) e al vettore (B).
void JFET::getStamps(
    int num_total_equations, double dt,
    const std::vector<double>& x,
    const std::vector<double>& prev_solution,
    double time,
    std::vector<std::vector<double>>& A,
    std::vector<double>& B
) {
    // Ottieni gli indici globali per i nodi del JFET.
    // Assumiamo che node_ids_ contenga [idx_D, idx_G, idx_S]
    // Questo è un cambiamento rispetto al tuo JFET.cpp originale che usava node_gate, node_drain, node_source come membri.
    // La classe base Component ora gestisce node_ids_ come vettore di int.
    // Assumiamo l'ordine dei nodi: Drain (0), Gate (1), Source (2)
    int idx_D = node_ids_[0];
    int idx_G = node_ids_[1];
    int idx_S = node_ids_[2];

    // Ottieni le tensioni correnti dei nodi dal vettore `x`
    // Il nodo 0 è il ground, quindi la sua tensione è 0.0
    double V_D_curr = (idx_D == 0) ? 0.0 : x[idx_D];
    double V_G_curr = (idx_G == 0) ? 0.0 : x[idx_G];
    double V_S_curr = (idx_S == 0) ? 0.0 : x[idx_S];

    // Calcola VGS e VDS correnti
    double VGS_curr = V_G_curr - V_S_curr;
    double VDS_curr = V_D_curr - V_S_curr; // Non usato direttamente nel modello JFET Level 1 per Id, ma utile per debug o modelli più complessi

    // Calcola gm e la corrente di drain equivalente basata su prev_VGS_
    double gm = calculateTransconductance(prev_VGS_);
    double Id_prev = calculateDrainCurrent(prev_VGS_);

    // Il modello linearizzato della corrente è: I_D = gm * VGS + I_eq
    // dove I_eq = Id_prev - gm * prev_VGS_
    // E VGS = V(idx_G) - V(idx_S)

    // Applica gli stamps del termine gm * VGS: gm * (V(idx_G) - V(idx_S))
    // Questa corrente fluisce da Drain (idx_D) a Source (idx_S) per N-canale.
    // È sottratta dall'equazione del nodo Drain e aggiunta all'equazione del nodo Source.
    // Per P-canale, la direzione della corrente è invertita.

    if (type == N_CHANNEL) {
        A[idx_D][idx_G] -= gm; // Contributo da V(idx_G) a I_D
        A[idx_D][idx_S] += gm; // Contributo da V(idx_S) a I_D

        A[idx_S][idx_G] += gm; // Contributo da V(idx_G) a I_S (che è -I_D)
        A[idx_S][idx_S] -= gm; // Contributo da V(idx_S) a I_S
    } else { // P_CHANNEL
        // Per PMOS, il gm è positivo, ma la corrente fluisce da Source a Drain.
        // Quindi gli stamps sono invertiti rispetto all'NMOS.
        A[idx_D][idx_G] += gm; // Contributo da V(idx_G) a I_D (che è negativo, quindi +gm)
        A[idx_D][idx_S] -= gm; // Contributo da V(idx_S) a I_D

        A[idx_S][idx_G] -= gm; // Contributo da V(idx_G) a I_S
        A[idx_S][idx_S] += gm; // Contributo da V(idx_S) a I_S
    }

    // Applica il termine della sorgente di corrente costante: I_eq = Id_prev - gm * prev_VGS_
    double I_eq = Id_prev - gm * prev_VGS_;

    // Questa corrente fluisce da Drain (idx_D) a Source (idx_S) per N-canale.
    // È aggiunta all'equazione del nodo Drain e sottratta dall'equazione del nodo Source.
    if (type == N_CHANNEL) {
        B[idx_D] -= I_eq; // Corrente esce dal Drain
        B[idx_S] += I_eq; // Corrente entra nel Source
    } else { // P_CHANNEL
        B[idx_D] += I_eq; // Corrente entra nel Drain
        B[idx_S] -= I_eq; // Corrente esce dal Source
    }
}

// Aggiorna la variabile di stato interna (prev_VGS_) per il prossimo passo temporale.
void JFET::updateState(const std::vector<double>& current_solution,
                         const std::vector<double>& prev_solution,
                         double dt) {
    // Ottieni gli indici globali per i nodi del JFET.
    // Assumiamo l'ordine dei nodi: Drain (0), Gate (1), Source (2)
    int idx_G = node_ids_[1];
    int idx_S = node_ids_[2];

    // Controlla se `current_solution` è abbastanza grande da contenere le tensioni dei nodi.
    if (idx_G >= current_solution.size() || idx_S >= current_solution.size()) {
        std::cerr << "Error: Dimensione del vettore soluzione non corrispondente o indice nodo non valido per JFET " << name_ << " in updateState." << std::endl;
        return;
    }

    // Ottieni le tensioni correnti dei nodi.
    // Il nodo 0 è il ground, quindi la sua tensione è 0.0
    double V_G_curr = (idx_G == 0) ? 0.0 : current_solution[idx_G];
    double V_S_curr = (idx_S == 0) ? 0.0 : current_solution[idx_S];

    // Aggiorna il VGS precedente per il prossimo passo temporale
    prev_VGS_ = V_G_curr - V_S_curr;
}
