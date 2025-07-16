// components/Cabinet_bass_reflex.cpp
#include "Cabinet_bass_reflex.h"
#include <iostream> // Per i messaggi di errore
#include <limits>   // Per numeric_limits
#include <cmath>    // Per M_PI

// Definisci PI se M_PI non è disponibile in cmath su alcuni compilatori
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Costruttore per il componente Cabinet_bass_reflex.
 *
 * Inizializza il modello di cassa bass-reflex con il suo nome, i nodi collegati,
 * la frequenza di accordo, il fattore di qualità e la resistenza di carico equivalente.
 * Calcola anche i valori equivalenti di L e C e inizializza le variabili di stato.
 *
 * @param name Il nome univoco del componente.
 * @param node_names_str Un vettore di stringhe contenente i nomi dei due nodi collegati.
 * @param Fb La frequenza di accordo del sistema bass-reflex (Hz). Deve essere > 0.
 * @param Qb Il fattore di qualità del sistema bass-reflex a Fb. Deve essere > 0.
 * @param R_load La resistenza elettrica equivalente per lo smorzamento/carico (Ohm). Deve essere > 0.
 */
Cabinet_bass_reflex::Cabinet_bass_reflex(const std::string& name, const std::vector<std::string>& node_names_str,
                                         double Fb, double Qb, double R_load)
    : Component(name, node_names_str), // Chiama il costruttore della classe base
      Fb_(Fb), Qb_(Qb), R_load_(R_load),
      L_eq_(0.0), C_eq_(0.0), // Saranno calcolati sotto
      prev_IL_(0.0), prev_VC_(0.0) // Inizializza le variabili di stato per l'analisi transitoria
{
    // Validazione dei parametri di input
    if (Fb_ <= 0.0) {
        std::cerr << "Attenzione: Cabinet_bass_reflex " << name_ << " ha una frequenza di accordo (Fb) non positiva. Impostato a 1.0Hz." << std::endl;
        this->Fb_ = 1.0;
    }
    if (Qb_ <= 0.0) {
        std::cerr << "Attenzione: Cabinet_bass_reflex " << name_ << " ha un fattore di qualità (Qb) non positivo. Impostato a 0.707." << std::endl;
        this->Qb_ = 0.707; // Valore tipico per risposta Butterworth
    }
    if (R_load_ <= 0.0) {
        std::cerr << "Attenzione: Cabinet_bass_reflex " << name_ << " ha una resistenza di carico (R_load) non positiva. Impostato a 8.0 Ohm." << std::endl;
        this->R_load_ = 8.0; // Impedenza tipica dell'altoparlante
    }

    // Calcola L e C equivalenti per un circuito RLC parallelo
    // Frequenza angolare: omega_b = 2 * PI * Fb
    // Per un RLC parallelo:
    // Q = R_load * sqrt(C_eq / L_eq)
    // omega_b = 1 / sqrt(L_eq * C_eq)
    //
    // Da queste, deriviamo:
    // C_eq = Qb / (R_load * omega_b)
    // L_eq = R_load / (Qb * omega_b)

    double omega_b = 2.0 * M_PI * this->Fb_;

    if (omega_b == 0.0) { // Evita la divisione per zero se Fb era in qualche modo 0 nonostante il controllo
        std::cerr << "Errore: omega_b calcolata è zero per Cabinet_bass_reflex " << name_ << ". Impossibile calcolare L_eq/C_eq." << std::endl;
        L_eq_ = std::numeric_limits<double>::infinity(); // Indica stato non valido
        C_eq_ = std::numeric_limits<double>::infinity();
    } else {
        C_eq_ = this->Qb_ / (this->R_load_ * omega_b);
        L_eq_ = this->R_load_ / (this->Qb_ * omega_b);
    }

    std::cout << "Cabinet_bass_reflex " << name_ << " inizializzato con Fb=" << Fb_ << "Hz, Qb=" << Qb_
              << ", R_load=" << R_load_ << " Ohm. L_eq equivalente=" << L_eq_ << " H, C_eq=" << C_eq_ << " F." << std::endl;
}

/**
 * @brief Applica gli "stamps" del circuito RLC parallelo equivalente alla matrice MNA (A) e al vettore (B).
 *
 * Questo metodo utilizza i modelli del compagno di Eulero all'indietro per l'induttore e il condensatore
 * per applicare i loro contributi in base allo stato del passo temporale precedente.
 *
 * @param num_total_equations Dimensione totale della matrice MNA.
 * @param dt Passo temporale. Cruciale per i modelli del compagno.
 * @param x Vettore della soluzione corrente (non usato per lo stamping, ma per updateState).
 * @param prev_solution Vettore della soluzione al passo temporale precedente (usato per i modelli del compagno).
 * @param time Tempo attuale della simulazione.
 * @param A Riferimento alla matrice MNA.
 * @param B Riferimento al vettore delle sorgenti (RHS).
 */
void Cabinet_bass_reflex::getStamps(
    int num_total_equations, double dt,
    const std::vector<double>& x,
    const std::vector<double>& prev_solution,
    double time,
    std::vector<std::vector<double>>& A,
    std::vector<double>& B
) {
    // Ottieni gli indici globali per i due nodi collegati.
    // Assumiamo che node_ids_ contenga [node1, node2]
    if (node_ids_.size() != 2) {
        std::cerr << "Errore: Cabinet_bass_reflex " << name_ << " si aspetta 2 ID di nodo, ma ne ha " << node_ids_.size() << std::endl;
        return;
    }
    int idx1 = node_ids_[0];
    int idx2 = node_ids_[1];

    // Controlla la validità degli indici
    if (idx1 < 0 || idx1 >= num_total_equations || idx2 < 0 || idx2 >= num_total_equations) {
        std::cerr << "Errore: Indice di nodo non valido per Cabinet_bass_reflex " << name_ << std::endl;
        return;
    }

    // Assicurati che dt sia positivo per l'analisi transitoria
    if (dt <= 0.0) {
        std::cerr << "Errore: Il passo temporale (dt) deve essere positivo per lo stamping transitorio di Cabinet_bass_reflex " << name_ << "." << std::endl;
        return;
    }

    // --- Applica lo stamp del Resistore parallelo (R_load) ---
    // Corrente attraverso il resistore: I_R = (V(idx1) - V(idx2)) / R_load
    double conductance_R = 1.0 / R_load_;
    if (idx1 != 0) A[idx1][idx1] += conductance_R;
    if (idx2 != 0) A[idx2][idx2] += conductance_R;
    if (idx1 != 0 && idx2 != 0) {
        A[idx1][idx2] -= conductance_R;
        A[idx2][idx1] -= conductance_R;
    }

    // --- Applica lo stamp dell'Induttore parallelo (L_eq) usando il modello del compagno di Eulero all'indietro ---
    // V_L(t) = L_eq * dI_L/dt  =>  V_L(t) = L_eq * (I_L(t) - I_L(t-dt)) / dt
    // Riorganizzando per I_L(t): I_L(t) = (dt / L_eq) * V_L(t) + I_L(t-dt)
    // Questo è equivalente a una conduttanza G_L = dt / L_eq in parallelo con una sorgente di corrente I_srcL = I_L(t-dt)
    if (L_eq_ == 0.0) { // Gestisci il corto circuito ideale se L_eq è zero
        std::cerr << "Attenzione: Cabinet_bass_reflex " << name_ << " ha L_eq=0. Trattato come corto circuito." << std::endl;
        // Effettivamente aggiungi una conduttanza molto grande o gestisci come un corto.
        // Per ora, evitiamo la divisione per zero e saltiamo lo stamping dell'induttore.
    } else {
        double conductance_L = dt / L_eq_;
        double I_srcL = prev_IL_; // Corrente dell'induttore dal passo temporale precedente

        if (idx1 != 0) A[idx1][idx1] += conductance_L;
        if (idx2 != 0) A[idx2][idx2] += conductance_L;
        if (idx1 != 0 && idx2 != 0) {
            A[idx1][idx2] -= conductance_L;
            A[idx2][idx1] -= conductance_L;
        }

        // La sorgente di corrente I_srcL fluisce da idx2 a idx1 (se V_L = V(idx1)-V(idx2))
        // Quindi, aggiungiamo I_srcL al nodo da cui esce (idx2) e sottraiamo dal nodo in cui entra (idx1).
        // No, la corrente I_srcL è una corrente che fluisce attraverso la conduttanza equivalente.
        // Se I_L(t) = G_L * V_L(t) + I_srcL, e V_L = V(idx1) - V(idx2),
        // allora la corrente entra in idx1 ed esce da idx2.
        // Quindi, I_srcL viene aggiunta a B[idx1] e sottratta da B[idx2].
        if (idx1 != 0) B[idx1] += I_srcL;
        if (idx2 != 0) B[idx2] -= I_srcL;
    }

    // --- Applica lo stamp del Condensatore parallelo (C_eq) usando il modello del compagno di Eulero all'indietro ---
    // I_C(t) = C_eq * dV_C/dt  =>  I_C(t) = C_eq * (V_C(t) - V_C(t-dt)) / dt
    // Questo è equivalente a una conduttanza G_C = C_eq / dt in parallelo con una sorgente di corrente I_srcC = (C_eq / dt) * V_C(t-dt)
    if (C_eq_ == 0.0) { // Gestisci il circuito aperto ideale se C_eq è zero
        // Questo significa effettivamente nessuna corrente del condensatore, quindi nessuno stamp necessario.
    } else {
        double conductance_C = C_eq_ / dt;
        // V_C_prev è la tensione del condensatore dal passo temporale precedente (V(idx1)_prev - V(idx2)_prev)
        double V_C_prev = prev_VC_;
        double I_srcC = conductance_C * V_C_prev;

        if (idx1 != 0) A[idx1][idx1] += conductance_C;
        if (idx2 != 0) A[idx2][idx2] += conductance_C;
        if (idx1 != 0 && idx2 != 0) {
            A[idx1][idx2] -= conductance_C;
            A[idx2][idx1] -= conductance_C;
        }

        // La sorgente di corrente I_srcC fluisce da idx2 a idx1.
        // Quindi, aggiungiamo I_srcC a B[idx1] e sottraiamo da B[idx2].
        if (idx1 != 0) B[idx1] += I_srcC;
        if (idx2 != 0) B[idx2] -= I_srcC;
    }
}

/**
 * @brief Aggiorna le variabili di stato interne (prev_IL, prev_VC) per il prossimo passo temporale.
 *
 * Questo metodo viene chiamato dopo che il sistema MNA è stato risolto per il passo temporale corrente.
 * Utilizza la tensione appena calcolata ai capi del componente (V(idx1)-V(idx2)) e la corrente dell'induttore precedente (prev_IL)
 * per aggiornare lo stato per la prossima iterazione/passo temporale.
 *
 * @param current_solution Il vettore della soluzione corrente (contiene le tensioni ai nodi).
 * @param prev_solution Il vettore della soluzione precedente (non usato direttamente qui per l'aggiornamento).
 * @param dt Il passo temporale.
 */
void Cabinet_bass_reflex::updateState(const std::vector<double>& current_solution,
                                      const std::vector<double>& prev_solution,
                                      double dt) {
    // Ottieni gli indici globali per i due nodi collegati.
    int idx1 = node_ids_[0];
    int idx2 = node_ids_[1];

    // Calcola la tensione corrente ai capi del componente.
    // Il nodo 0 è il ground.
    double V_curr = ((idx1 == 0) ? 0.0 : current_solution[idx1]) -
                    ((idx2 == 0) ? 0.0 : current_solution[idx2]);

    // Aggiorna la tensione del condensatore per il prossimo passo temporale
    prev_VC_ = V_curr;

    // Aggiorna la corrente dell'induttore per il prossimo passo temporale usando la formula di Eulero all'indietro
    // I_L(t) = I_L(t-dt) + (dt / L_eq) * V_L(t)
    // Qui, V_curr è V_L(t)
    if (L_eq_ != 0.0) { // Evita la divisione per zero
        prev_IL_ += (dt / L_eq_) * V_curr;
    } else {
        // Se L_eq è 0, è un corto circuito, la corrente può essere infinita o indefinita.
        // Per la simulazione pratica, questo caso dovrebbe essere gestito dall'utente fornendo L_eq non zero.
        prev_IL_ = std::numeric_limits<double>::quiet_NaN(); // Indica stato non valido
    }
}
