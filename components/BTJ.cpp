// components/BTJ.cpp
#include "BTJ.h"
#include <iostream> // Per messaggi di errore
#include <limits>   // Per std::numeric_limits (es. quiet_NaN)

/**
 * @brief Costruttore per il componente BTJ.
 *
 * Inizializza un BJT NPN con il suo nome, i tre nodi collegati (Base, Collettore, Emettitore),
 * e i parametri chiave.
 *
 * @param name Il nome univoco del BJT.
 * @param node_base Il nome del nodo di Base.
 * @param node_collector Il nome del nodo di Collettore.
 * @param node_emitter Il nome del nodo di Emettitore.
 * @param Is La corrente di saturazione (Ampere). Deve essere > 0.
 * @param beta_F Il guadagno di corrente diretto in configurazione a emettitore comune. Deve essere > 0.
 * @param Vt La tensione termica (Volt). Il valore predefinito è 0.02585V (circa a 300K). Deve essere > 0.
 * @param tolerance_percent Percentuale di tolleranza opzionale.
 */
BTJ::BTJ(const std::string& name,
         const std::string& node_base, const std::string& node_collector, const std::string& node_emitter,
         double Is, double beta_F, double Vt, double tolerance_percent)
    // Il costruttore della classe base Component accetta node1 e node2.
    // Per un dispositivo a 3 terminali, possiamo mappare node1 alla Base e node2 al Collettore,
    // e gestire node_emitter separatamente. I node1 e node2 della classe base
    // non sono strettamente usati per le connessioni interne del BJT, ma per la sua
    // presenza complessiva nel circuito.
    : Component(name, node_base, node_collector, tolerance_percent), // node1=base, node2=collector
      node_base(node_base), node_collector(node_collector), node_emitter(node_emitter),
      Is(Is), Vt(Vt), beta_F(beta_F), alpha_F(0.0), // alpha_F calcolato sotto
      prev_VBE(0.0), prev_VBC(0.0) // Inizializza le variabili di stato
{
    // Un modello di BJT tipicamente non richiede variabili ausiliarie aggiuntive oltre
    // alle tensioni dei nodi quando si utilizzano modelli equivalenti per le non linearità.
    // La sorgente di corrente dipendente viene gestita direttamente tramite stamping nella matrice A.
    setNumAuxiliaryVariables(0);

    // Valida i parametri di input
    if (Is <= 0.0) {
        std::cerr << "Attenzione: BTJ " << name << " ha una corrente di saturazione (Is) non positiva. Impostazione a 1e-15A." << std::endl;
        this->Is = 1e-15; // Valore tipico predefinito
    }
    if (beta_F <= 0.0) {
        std::cerr << "Attenzione: BTJ " << name << " ha un beta diretto (beta_F) non positivo. Impostazione a 100." << std::endl;
        this->beta_F = 100.0;
    }
    if (Vt <= 0.0) {
        std::cerr << "Attenzione: BTJ " << name << " ha una tensione termica (Vt) non positiva. Impostazione a 0.02585V." << std::endl;
        this->Vt = 0.02585; // Circa a 300K
    }

    // Calcola alpha_F da beta_F
    alpha_F = beta_F / (1.0 + beta_F);
    if (beta_F == std::numeric_limits<double>::infinity()) { // Gestisce beta infinito
        alpha_F = 1.0;
    }

    std::cout << "BTJ " << name << " inizializzato con Is=" << Is << "A, beta_F=" << beta_F
              << ", Vt=" << Vt << "V. Alpha_F=" << alpha_F << "." << std::endl;
}

/**
 * @brief Funzione helper per calcolare i parametri del modello equivalente del diodo.
 *
 * Calcola la conduttanza equivalente (G_eq) e la sorgente di corrente (I_src)
 * per un diodo utilizzando un modello linearizzato basato sulla tensione precedente.
 *
 * Corrente del diodo: $I_D = I_S \cdot (e^{V_D / V_T} - 1)$
 * Modello linearizzato: $I_D = G_{eq} \cdot V_D - I_{src}$
 * Dove:
 * $G_{eq} = (dI_D / dV_D)$ valutato a $V_{D\_prev} = (I_S / V_T) \cdot e^{V_{D\_prev} / V_T}$
 * $I_{src} = G_{eq} \cdot V_{D\_prev} - I_S \cdot (e^{V_{D\_prev} / V_T} - 1)$
 *
 * @param V_prev La tensione attraverso il diodo dal passo temporale precedente.
 * @param G_eq Output: La conduttanza equivalente del diodo.
 * @param I_src Output: La sorgente di corrente equivalente del diodo.
 */
void BTJ::calculateDiodeCompanionModel(double V_prev, double& G_eq, double& I_src) {
    // Evita esponenti molto grandi per prevenire l'overflow
    double exp_term;
    if (V_prev / Vt > 700.0) { // Valore arbitrario grande per prevenire overflow per exp()
        exp_term = std::numeric_limits<double>::max();
    } else if (V_prev / Vt < -700.0) { // Valore arbitrario piccolo per underflow
        exp_term = 0.0;
    } else {
        exp_term = std::exp(V_prev / Vt);
    }

    G_eq = (Is / Vt) * exp_term;
    I_src = G_eq * V_prev - Is * (exp_term - 1.0);

    // Limita G_eq a un massimo ragionevole per prevenire instabilità numerica
    // Una conduttanza molto grande può causare problemi nell'inversione della matrice.
    const double MAX_G_EQ = 1e12; // Esempio: 1e12 Siemens (equivalente a 1e-12 Ohm)
    if (G_eq > MAX_G_EQ) {
        G_eq = MAX_G_EQ;
    }
}

/**
 * @brief Applica gli "stamps" del modello Ebers-Moll del BJT alla matrice MNA (A) e al vettore (b).
 *
 * Questo metodo utilizza modelli equivalenti linearizzati per i diodi Base-Emettitore e Base-Collettore
 * basati sulle tensioni del passo temporale precedente. Applica anche lo stamp della sorgente
 * di corrente dipendente (alpha_F * I_F) dal Collettore all'Emettitore.
 *
 * @param A La matrice MNA a cui vengono applicati gli stamps.
 * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
 * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo (usata per recuperare le tensioni della soluzione precedente).
 * @param prev_solution La soluzione dal passo temporale precedente (usata per i modelli equivalenti).
 * @param time Il tempo di simulazione corrente (non usato direttamente per lo stamping).
 * @param dt La dimensione del passo temporale (non usato direttamente per questo modello statico, ma influenza implicitamente la soluzione precedente).
 */
void BTJ::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Ottieni gli indici globali per i nodi del BJT.
    // Si assume che getNodeIndex sia disponibile e mappi correttamente i nomi dei nodi agli indici della matrice.
    int idx_B = getNodeIndex(node_base);
    int idx_C = getNodeIndex(node_collector);
    int idx_E = getNodeIndex(node_emitter);

    // Controlla la validità degli indici
    if (idx_B == -1 || idx_C == -1 || idx_E == -1) {
        std::cerr << "Errore: Indice di nodo non valido per BTJ " << name << std::endl;
        return;
    }

    double G_BE, I_src_BE;
    double G_BC, I_src_BC;

    // --- 1. Modello Equivalente del Diodo Base-Emettitore ---
    // prev_VBE è la tensione attraverso la giunzione BE dal passo temporale precedente.
    calculateDiodeCompanionModel(prev_VBE, G_BE, I_src_BE);

    // Applica G_BE tra Base (idx_B) ed Emettitore (idx_E)
    A(idx_B, idx_B) += G_BE;
    A(idx_E, idx_E) += G_BE;
    A(idx_B, idx_E) -= G_BE;
    A(idx_E, idx_B) -= G_BE;

    // Applica I_src_BE come sorgente di corrente dall'Emettitore alla Base
    b(idx_B) += I_src_BE;
    b(idx_E) -= I_src_BE;

    // --- 2. Modello Equivalente del Diodo Base-Collettore ---
    // prev_VBC è la tensione attraverso la giunzione BC dal passo temporale precedente.
    calculateDiodeCompanionModel(prev_VBC, G_BC, I_src_BC);

    // Applica G_BC tra Base (idx_B) e Collettore (idx_C)
    A(idx_B, idx_B) += G_BC;
    A(idx_C, idx_C) += G_BC;
    A(idx_B, idx_C) -= G_BC;
    A(idx_C, idx_B) -= G_BC;

    // Applica I_src_BC come sorgente di corrente dal Collettore alla Base
    b(idx_B) += I_src_BC;
    b(idx_C) -= I_src_BC;

    // --- 3. Sorgente di Corrente Dipendente (alpha_F * I_F) dal Collettore all'Emettitore ---
    // I_F è la corrente diretta attraverso il diodo BE.
    // Dal modello equivalente del diodo BE: $I_F = G_{BE} \cdot V_{BE} - I_{src\_BE}$
    // Dove $V_{BE} = V(idx\_B) - V(idx\_E)$
    // Quindi la sorgente di corrente dipendente è $I_{dep} = \alpha_F \cdot (G_{BE} \cdot (V(idx\_B) - V(idx\_E)) - I_{src\_BE})$

    // Applica la parte controllata dalla tensione ($\alpha_F \cdot G_{BE} \cdot (V(idx\_B) - V(idx\_E)))$
    // Questa corrente fluisce dal Collettore (idx_C) all'Emettitore (idx_E).
    // Sottrae dall'equazione del nodo Collettore e aggiunge all'equazione del nodo Emettitore.
    A(idx_C, idx_B) -= alpha_F * G_BE; // Contributo da V(idx_B) a I_C
    A(idx_C, idx_E) += alpha_F * G_BE; // Contributo da V(idx_E) a I_C

    A(idx_E, idx_B) += alpha_F * G_BE; // Contributo da V(idx_B) a I_E
    A(idx_E, idx_E) -= alpha_F * G_BE; // Contributo da V(idx_E) a I_E

    // Applica la parte costante ($\alpha_F \cdot I_{src\_BE}$)
    // Questa corrente fluisce dal Collettore (idx_C) all'Emettitore (idx_E).
    b(idx_C) += alpha_F * I_src_BE; // Aggiunge all'equazione del nodo Collettore
    b(idx_E) -= alpha_F * I_src_BE; // Sottrae dall'equazione del nodo Emettitore
}

/**
 * @brief Aggiorna le variabili di stato interne (prev_VBE, prev_VBC) per il passo temporale successivo.
 *
 * Questo metodo viene chiamato dopo che il sistema MNA è stato risolto per il passo temporale corrente.
 * Estrae le attuali tensioni Base-Emettitore e Base-Collettore dal vettore della soluzione
 * (x_current_guess) e le memorizza per l'uso nei modelli equivalenti del passo temporale successivo.
 *
 * @param v_curr La tensione corrente attraverso il componente (non usata direttamente per gli stati interni dei diodi).
 * @param i_curr La corrente che attraversa il componente (non usata direttamente).
 */
void BTJ::updateState(double v_curr, double i_curr) {
    // I parametri `v_curr` e `i_curr` qui sono per il componente nel suo complesso,
    // il che non è direttamente significativo per un dispositivo a 3 terminali come un BJT.
    // Invece, dobbiamo accedere alle tensioni dei nodi risolte dal vettore completo della soluzione.
    // Questo tipicamente avviene nel ciclo principale del Circuit dopo la risoluzione di A*x=b.
    // Si assume che x_current_guess (che è il vettore della soluzione `x` dal passo temporale corrente)
    // sia disponibile e contenga le ultime tensioni dei nodi.

    // Ottieni gli indici globali per i nodi del BJT.
    int idx_B = getNodeIndex(node_base);
    int idx_C = getNodeIndex(node_collector);
    int idx_E = getNodeIndex(node_emitter);

    // Recupera le tensioni attuali dei nodi dal vettore della soluzione.
    // Questo presuppone che la classe Circuit passi il *vettore completo della soluzione*
    // come `x_current_guess` a `updateState` dopo `getStamps` e la risoluzione.
    // Se `updateState` viene chiamato con `v_curr` come tensione tra `node1` e `node2`
    // della classe base Component, allora questo deve essere rivalutato.
    // Assumendo che `x_current_guess` sia il vettore completo della soluzione `x`.

    // Controlla se `x_current_guess` è abbastanza grande da contenere le tensioni dei nodi.
    // Il numero di tensioni dei nodi è tipicamente `num_nodes`.
    // Questo controllo è una salvaguardia; l'integrazione corretta con il risolutore del Circuito è fondamentale.
    if (idx_B < 0 || idx_B >= x_current_guess.size() ||
        idx_C < 0 || idx_C >= x_current_guess.size() ||
        idx_E < 0 || idx_E >= x_current_guess.size()) {
        std::cerr << "Errore: Dimensione del vettore della soluzione non corrispondente o indice di nodo non valido per BTJ " << name << " in updateState." << std::endl;
        return;
    }

    double V_B_curr = x_current_guess(idx_B);
    double V_C_curr = x_current_guess(idx_C);
    double V_E_curr = x_current_guess(idx_E);

    // Aggiorna le tensioni delle giunzioni precedenti per il passo temporale successivo
    prev_VBE = V_B_curr - V_E_curr;
    prev_VBC = V_B_curr - V_C_curr;
}
