// components/BTJ.h
#ifndef BTJ_H
#define BTJ_H

#include "Component.h" // Include la classe base Component
#include <Eigen/Dense> // Per le operazioni con matrici e vettori Eigen
#include <string>      // Per std::string
#include <cmath>       // Per std::exp

/**
 * @class BTJ
 * @brief Rappresenta un Transistor a Giunzione Bipolare (BJT) utilizzando un modello semplificato di Ebers-Moll.
 *
 * Questa classe modella un BJT NPN utilizzando un modello semplificato di Ebers-Moll, che consiste in
 * due diodi (Base-Emettitore e Base-Collettore) e una sorgente di corrente controllata in corrente (CCCS)
 * dal Collettore all'Emettitore.
 *
 * Per l'analisi transitoria, le caratteristiche non lineari dei diodi sono linearizzate utilizzando
 * modelli equivalenti (companion models) basati sulla tensione dal passo temporale precedente. Questo approccio
 * è adatto per una simulazione di base nel dominio del tempo ma potrebbe richiedere passi temporali piccoli
 * o un risolutore iterativo esterno per non linearità forti.
 *
 * I parametri del modello includono la corrente di saturazione (Is), la tensione termica (Vt) e
 * il guadagno di corrente diretto (beta_F). Il guadagno di corrente inverso (beta_R) è trascurato
 * per semplicità in questo modello.
 */
class BTJ : public Component {
public:
    // Nodi specifici del BJT
    std::string node_base;    ///< Nome del nodo di Base.
    std::string node_collector; ///< Nome del nodo di Collettore.
    std::string node_emitter; ///< Nome del nodo di Emettitore.

    // Parametri del BJT
    double Is;       ///< Corrente di saturazione (Ampere).
    double Vt;       ///< Tensione termica (Volt), tipicamente ~25.85mV a 300K.
    double beta_F;   ///< Guadagno di corrente in configurazione a emettitore comune (h_FE) diretto.
    double alpha_F;  ///< Guadagno di corrente in configurazione a base comune (derivato da beta_F).

    // Variabili di stato per i modelli equivalenti (tensioni attraverso le giunzioni dal passo temporale precedente)
    double prev_VBE; ///< Tensione Base-Emettitore dal passo temporale precedente.
    double prev_VBC; ///< Tensione Base-Collettore dal passo temporale precedente.

    /**
     * @brief Costruttore per il componente BTJ.
     *
     * Inizializza un BJT NPN con il suo nome, i tre nodi collegati (Base, Collettore, Emettitore)
     * e i parametri chiave.
     *
     * @param name Il nome univoco del BJT.
     * @param node_base Il nome del nodo di Base.
     * @param node_collector Il nome del nodo di Collettore.
     * @param node_emitter Il nome del nodo di Emettitore.
     * @param Is La corrente di saturazione (Ampere). Deve essere > 0.
     * @param beta_F Il guadagno di corrente diretto in configurazione a emettitore comune. Deve essere > 0.
     * @param Vt La tensione termica (Volt). Il valore predefinito è 0.02585V (circa a 300K). Deve essere > 0.
     * @param tolerance_percent Percentuale di tolleranza opzionale (non usata direttamente per questo modello).
     */
    BTJ(const std::string& name,
        const std::string& node_base, const std::string& node_collector, const std::string& node_emitter,
        double Is, double beta_F, double Vt = 0.02585, double tolerance_percent = 0.0);

    /**
     * @brief Crea una copia profonda dell'oggetto BTJ.
     * @return Un puntatore a un nuovo oggetto BTJ, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new BTJ(*this); }

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
     * @param dt La dimensione del passo temporale (non usata direttamente per questo modello statico, ma influenza implicitamente la soluzione precedente).
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Aggiorna le variabili di stato interne (prev_VBE, prev_VBC) per il passo temporale successivo.
     *
     * Questo metodo viene chiamato dopo che il sistema MNA è stato risolto per il passo temporale corrente.
     * Estrae le attuali tensioni Base-Emettitore e Base-Collettore dal vettore della soluzione
     * e le memorizza per l'uso nei modelli equivalenti del passo temporale successivo.
     *
     * @param v_curr La tensione corrente attraverso il componente (non usata direttamente per gli stati interni dei diodi).
     * @param i_curr La corrente che attraversa il componente (non usata direttamente).
     */
    void updateState(double v_curr, double i_curr) override;

private:
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
    void calculateDiodeCompanionModel(double V_prev, double& G_eq, double& I_src);
};

#endif // BTJ_H
