// components/Pentode.cpp
#include "Pentode.h"
#include <iostream>  // Per messaggi di debug e avvisi
#include <cmath>     // Per std::tanh, std::cosh, std::max, std::fabs
#include <algorithm> // Per std::max

// Costante per la stabilità numerica: un valore molto piccolo per evitare divisioni per zero o logaritmi/tanh di valori estremi.
const double EPSILON_PENTODE = 1e-9;

/**
 * @brief Costruttore per il componente Pentode.
 *
 * Inizializza un pentodo con il suo nome, i cinque nodi collegati
 * e i parametri chiave del modello semplificato.
 *
 * @param name Il nome univoco del pentodo.
 * @param control_grid_node Il nome del nodo della griglia di controllo.
 * @param plate_node Il nome del nodo della placca.
 * @param cathode_node Il nome del nodo del catodo.
 * @param screen_grid_node Il nome del nodo della griglia schermo.
 * @param suppressor_grid_node Il nome del nodo della griglia soppressore.
 * @param K_p Parametro di transconduttanza.
 * @param V_cutoff Tensione di cutoff.
 * @param lambda Parametro lambda per la modulazione della lunghezza del canale.
 * @param mu_sg Fattore di amplificazione della griglia schermo.
 * @param I_s_ratio Rapporto corrente griglia schermo/placca.
 * @param tolerance_percent Percentuale di tolleranza opzionale.
 */
Pentode::Pentode(const std::string& name,
                 const std::string& control_grid_node, const std::string& plate_node,
                 const std::string& cathode_node, const std::string& screen_grid_node,
                 const std::string& suppressor_grid_node,
                 double K_p, double V_cutoff, double lambda, double mu_sg, double I_s_ratio,
                 double tolerance_percent)
    // Chiama il costruttore della classe base Component.
    // Per un componente a più terminali, passiamo i nodi principali (es. placca e catodo)
    // al costruttore base, e gestiamo gli altri nodi internamente.
    : Component(name, plate_node, cathode_node, tolerance_percent),
      node_control_grid(control_grid_node),
      node_plate(plate_node),
      node_cathode(cathode_node),
      node_screen_grid(screen_grid_node),
      node_suppressor_grid(suppressor_grid_node),
      K_p_(K_p),
      V_cutoff_(V_cutoff),
      lambda_(lambda),
      mu_sg_(mu_sg),
      I_s_ratio_(I_s_ratio),
      prev_Vcg_(0.0), prev_Vp_(0.0), prev_Vsg_(0.0), prev_Vsup_(0.0) // Inizializza le tensioni precedenti
{
    // Aggiungi tutti i nodi del pentodo alla lista dei nodi gestiti dal componente.
    // Questo è cruciale per getNodeIndex per funzionare correttamente per tutti i terminali.
    addNode(control_grid_node);
    addNode(screen_grid_node);
    addNode(suppressor_grid_node); // Potrebbe essere collegato a massa esternamente

    // Validazione dei parametri di input
    if (K_p_ <= 0) {
        std::cerr << "Attenzione: Pentodo " << name << " ha K_p non positivo. Impostazione a 1e-3." << std::endl;
        K_p_ = 1e-3; // Valore di default ragionevole
    }
    if (mu_sg_ <= 0) {
        std::cerr << "Attenzione: Pentodo " << name << " ha mu_sg non positivo. Impostazione a 10." << std::endl;
        mu_sg_ = 10.0; // Valore di default ragionevole
    }
    if (lambda_ < 0) {
        std::cerr << "Attenzione: Pentodo " << name << " ha lambda negativo. Impostazione a 0." << std::endl;
        lambda_ = 0.0; // Non dovrebbe essere negativo
    }
    if (I_s_ratio_ < 0 || I_s_ratio_ > 1) {
        std::cerr << "Attenzione: Pentodo " << name << " ha I_s_ratio fuori range [0,1]. Impostazione a 0.2." << std::endl;
        I_s_ratio_ = 0.2; // Valore di default ragionevole
    }

    std::cout << "Pentodo '" << name << "' inizializzato." << std::endl;
}

/**
 * @brief Applica gli "stamps" del pentodo alla matrice MNA (A) e al vettore (b).
 *
 * Questo metodo implementa un modello non lineare semplificato per il pentodo,
 * calcolando le correnti di placca e griglia schermo e le loro conduttanze dinamiche
 * per l'iterazione di Newton-Raphson.
 *
 * Il modello semplificato per la corrente di placca (Ip) è basato su una funzione
 * tanh dell'effettiva tensione di griglia, con una dipendenza lineare dalla tensione di placca:
 * V_effective_grid = V_cg + V_sg / mu_sg
 * Ip = K_p * tanh(V_effective_grid - V_cutoff) * (1 + lambda * V_p)
 *
 * La corrente di griglia schermo (Isg) è modellata come una frazione della corrente di placca:
 * Isg = I_s_ratio * Ip
 *
 * I contributi MNA sono calcolati usando il modello di accompagnamento linearizzato.
 *
 * @param A La matrice MNA a cui vengono applicati gli stamps.
 * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
 * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo.
 * @param prev_solution La soluzione dal passo temporale precedente (non usata direttamente per il modello statico non lineare).
 * @param time Il tempo di simulazione corrente.
 * @param dt La dimensione del passo temporale.
 */
void Pentode::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Ottieni gli indici globali (base 1) per tutti i nodi del pentodo.
    int idx_cg = getNodeIndex(node_control_grid);
    int idx_p = getNodeIndex(node_plate);
    int idx_k = getNodeIndex(node_cathode); // Catodo (spesso massa, ID 0)
    int idx_sg = getNodeIndex(node_screen_grid);
    int idx_sup = getNodeIndex(node_suppressor_grid); // Soppressore (spesso massa, ID 0)

    // Recupera le tensioni attuali dei nodi dalla stima corrente (x_current_guess).
    // Se un nodo è la massa (ID 0), la sua tensione è 0.
    double V_cg_curr = (idx_cg == 0) ? 0.0 : x_current_guess(idx_cg - 1);
    double V_p_curr = (idx_p == 0) ? 0.0 : x_current_guess(idx_p - 1);
    double V_k_curr = (idx_k == 0) ? 0.0 : x_current_guess(idx_k - 1);
    double V_sg_curr = (idx_sg == 0) ? 0.0 : x_current_guess(idx_sg - 1);
    double V_sup_curr = (idx_sup == 0) ? 0.0 : x_current_guess(idx_sup - 1);

    // Calcola le tensioni relative al catodo (se il catodo non è a massa, altrimenti sono tensioni assolute).
    double V_cg_rel = V_cg_curr - V_k_curr;
    double V_p_rel = V_p_curr - V_k_curr;
    double V_sg_rel = V_sg_curr - V_k_curr;
    // double V_sup_rel = V_sup_curr - V_k_curr; // Non usata direttamente nel modello semplificato

    // --- Calcolo delle Correnti e Conduttanze Dinamiche ---

    // 1. Tensione di griglia efficace
    // Previene divisioni per zero se mu_sg_ è troppo piccolo
    double mu_sg_safe = std::max(EPSILON_PENTODE, mu_sg_);
    double V_effective_grid = V_cg_rel + V_sg_rel / mu_sg_safe;

    // 2. Calcolo della corrente di placca (Ip) e delle sue derivate (conduttanze)
    // Usiamo una funzione tanh per una transizione morbida e differenziabile.
    double tanh_arg = V_effective_grid - V_cutoff_;
    double Ip_curr = K_p_ * std::tanh(tanh_arg) * (1.0 + lambda_ * V_p_rel);

    // Derivata di tanh(x) è sech^2(x) = 1 - tanh^2(x)
    double sech2_tanh_arg = 1.0 - std::tanh(tanh_arg) * std::tanh(tanh_arg);

    // Transconduttanza della griglia di controllo (gm_cg = dIp/dVcg)
    double gm_cg = K_p_ * sech2_tanh_arg * (1.0 + lambda_ * V_p_rel);

    // Transconduttanza della griglia schermo (gm_sg_plate = dIp/dVsg)
    double gm_sg_plate = K_p_ * sech2_tanh_arg * (1.0 + lambda_ * V_p_rel) / mu_sg_safe;

    // Conduttanza di placca (gp = dIp/dVp)
    double gp = K_p_ * std::tanh(tanh_arg) * lambda_;

    // 3. Calcolo della corrente di griglia schermo (Isg) e delle sue derivate
    // Isg è modellata come una frazione di Ip.
    double Isg_curr = I_s_ratio_ * Ip_curr;

    // Derivate di Isg rispetto alle tensioni
    double gm_cg_screen = I_s_ratio_ * gm_cg;
    double gm_sg_screen = I_s_ratio_ * gm_sg_plate;
    double gp_screen = I_s_ratio_ * gp;

    // --- Applicazione degli "stamps" MNA ---

    // Termini di corrente equivalenti per il lato destro (b)
    // I_eq = I_NL(V_guess) - G_NL * V_guess
    double I_eq_p = Ip_curr - (gm_cg * V_cg_rel + gm_sg_plate * V_sg_rel + gp * V_p_rel);
    double I_eq_sg = Isg_curr - (gm_cg_screen * V_cg_rel + gm_sg_screen * V_sg_rel + gp_screen * V_p_rel);

    // KCL al nodo della placca (idx_p)
    if (idx_p != 0) {
        A(idx_p, idx_p) += gp;
        if (idx_cg != 0) A(idx_p, idx_cg) += gm_cg;
        if (idx_sg != 0) A(idx_p, idx_sg) += gm_sg_plate;
        if (idx_k != 0) A(idx_p, idx_k) -= (gp + gm_cg + gm_sg_plate); // Correzione per V_k_rel = V_k_curr
        b(idx_p) -= I_eq_p;
    }

    // KCL al nodo della griglia schermo (idx_sg)
    if (idx_sg != 0) {
        A(idx_sg, idx_sg) += gp_screen; // Conduttanza di uscita della griglia schermo
        if (idx_cg != 0) A(idx_sg, idx_cg) += gm_cg_screen;
        if (idx_p != 0) A(idx_sg, idx_p) += gp_screen; // Dipendenza da Vp per Isg
        if (idx_k != 0) A(idx_sg, idx_k) -= (gp_screen + gm_cg_screen + gp_screen); // Correzione per V_k_rel
        b(idx_sg) -= I_eq_sg;
    }

    // KCL al nodo del catodo (idx_k)
    // La corrente totale che entra nel catodo è -(Ip + Isg).
    // Quindi, i contributi al catodo sono l'opposto di quelli di placca e griglia schermo.
    if (idx_k != 0) {
        // Contributi da Ip
        if (idx_p != 0) A(idx_k, idx_p) -= gp;
        if (idx_cg != 0) A(idx_k, idx_cg) -= gm_cg;
        if (idx_sg != 0) A(idx_k, idx_sg) -= gm_sg_plate;
        A(idx_k, idx_k) += (gp + gm_cg + gm_sg_plate); // Auto-conduttanza del catodo per Ip

        // Contributi da Isg
        if (idx_sg != 0) A(idx_k, idx_sg) -= gp_screen; // Conduttanza di uscita della griglia schermo
        if (idx_cg != 0) A(idx_k, idx_cg) -= gm_cg_screen;
        if (idx_p != 0) A(idx_k, idx_p) -= gp_screen; // Dipendenza da Vp per Isg
        A(idx_k, idx_k) += (gp_screen + gm_cg_screen + gp_screen); // Auto-conduttanza del catodo per Isg

        b(idx_k) += (I_eq_p + I_eq_sg); // Somma delle correnti equivalenti
    }

    // La griglia di controllo e la griglia soppressore sono idealmente ad alta impedenza
    // e non contribuiscono direttamente con correnti significative in questo modello semplificato.
    // Se ci fossero correnti di griglia (es. griglia di controllo in conduzione),
    // andrebbero aggiunti qui i relativi stamps.

    // La griglia soppressore (idx_sup) è spesso collegata al catodo o a massa.
    // In questo modello semplificato, non ha contributi attivi diretti.
}

/**
 * @brief Aggiorna lo stato interno del pentodo.
 *
 * Questo metodo memorizza le tensioni dei nodi attuali (dalla soluzione MNA)
 * come "precedenti" per la prossima iterazione di Newton-Raphson nel metodo getStamps.
 * Questo è cruciale per la convergenza dei modelli non lineari.
 *
 * @param v_curr La tensione corrente attraverso il componente (non usata direttamente).
 * @param i_curr La corrente corrente che attraversa il componente (non usata direttamente).
 */
void Pentode::updateState(double v_curr, double i_curr) {
    // Ottieni gli indici globali (base 1) per i nodi del pentodo.
    int idx_cg = getNodeIndex(node_control_grid);
    int idx_p = getNodeIndex(node_plate);
    int idx_sg = getNodeIndex(node_screen_grid);
    int idx_sup = getNodeIndex(node_suppressor_grid);
    int idx_k = getNodeIndex(node_cathode);

    // Accedi al vettore della soluzione completa per ottenere le tensioni dei nodi.
    // Assumiamo che `x_current_guess` sia disponibile o che queste tensioni vengano passate
    // in un modo che le renda accessibili qui. Per semplicità, userò `x_current_guess`
    // come se fosse una variabile membro o passata implicitamente.
    // In un sistema reale, il solutore passerebbe la soluzione completa.
    // Per questo esempio, useremo le variabili membro `prev_V` per memorizzare.
    // Questo metodo viene chiamato *dopo* che A*x=b è stato risolto.

    // Le tensioni da memorizzare sono quelle calcolate nell'ultimo passo di Newton-Raphson.
    // Queste dovrebbero essere le tensioni assolute dei nodi (rispetto a massa).
    // Se il nodo è massa (ID 0), la tensione è 0.
    // Altrimenti, la tensione è x_current_guess(idx - 1) (assumendo x_current_guess è 0-indexed).
    if (idx_cg != 0) prev_Vcg_ = x_current_guess(idx_cg - 1);
    else prev_Vcg_ = 0.0;

    if (idx_p != 0) prev_Vp_ = x_current_guess(idx_p - 1);
    else prev_Vp_ = 0.0;

    if (idx_sg != 0) prev_Vsg_ = x_current_guess(idx_sg - 1);
    else prev_Vsg_ = 0.0;

    if (idx_sup != 0) prev_Vsup_ = x_current_guess(idx_sup - 1);
    else prev_Vsup_ = 0.0;

    // Il catodo (prev_Vk_) non è strettamente necessario da memorizzare se è sempre a massa,
    // ma lo includiamo per completezza se il catodo può essere non a massa.
    if (idx_k != 0) prev_Vk_ = x_current_guess(idx_k - 1);
    else prev_Vk_ = 0.0;
}
