// components/LDR.cpp
#include "LDR.h"
#include <iostream> // Per messaggi di errore
#include <limits>   // Per std::numeric_limits

/**
 * @brief Costruttore per il componente LDR.
 *
 * Inizializza un LDR con il suo nome, i due nodi collegati e i parametri del modello.
 *
 * @param name Il nome univoco dell'LDR.
 * @param node1 Il nome del primo nodo.
 * @param node2 Il nome del secondo nodo.
 * @param R0 La resistenza di riferimento (Ohm) a L0 lux. Deve essere > 0.
 * @param L0 L'intensità luminosa di riferimento (lux). Deve essere > 0.
 * @param gamma Il coefficiente di sensibilità alla luce. Deve essere > 0.
 * @param initial_light_intensity L'intensità luminosa iniziale (lux). Deve essere >= 0.
 * @param tolerance_percent Percentuale di tolleranza opzionale.
 */
LDR::LDR(const std::string& name,
         const std::string& node1, const std::string& node2,
         double R0, double L0, double gamma, double initial_light_intensity,
         double tolerance_percent)
    : Component(name, node1, node2, tolerance_percent),
      R0(R0), L0(L0), gamma(gamma), current_light_intensity(initial_light_intensity)
{
    // Un resistore non richiede variabili ausiliarie aggiuntive.
    setNumAuxiliaryVariables(0);

    // Validazione dei parametri di input
    if (R0 <= 0.0) {
        std::cerr << "Attenzione: LDR " << name << " ha una resistenza di riferimento (R0) non positiva. Impostazione a 1000 Ohm." << std::endl;
        this->R0 = 1000.0;
    }
    if (L0 <= 0.0) {
        std::cerr << "Attenzione: LDR " << name << " ha un'intensità luminosa di riferimento (L0) non positiva. Impostazione a 10 lux." << std::endl;
        this->L0 = 10.0;
    }
    if (gamma <= 0.0) {
        std::cerr << "Attenzione: LDR " << name << " ha un coefficiente gamma non positivo. Impostazione a 0.7." << std::endl;
        this->gamma = 0.7;
    }
    if (initial_light_intensity < 0.0) {
        std::cerr << "Attenzione: LDR " << name << " ha un'intensità luminosa iniziale negativa. Impostazione a 0 lux." << std::endl;
        this->current_light_intensity = 0.0;
    }

    std::cout << "LDR " << name << " inizializzato con R0=" << R0 << " Ohm, L0=" << L0
              << " lux, gamma=" << gamma << ", luce iniziale=" << initial_light_intensity << " lux." << std::endl;
}

/**
 * @brief Applica gli "stamps" del modello LDR alla matrice MNA (A) e al vettore (b).
 *
 * Calcola la resistenza attuale dell'LDR in base a `current_light_intensity`
 * e applica lo stamp di un resistore.
 *
 * @param A La matrice MNA a cui vengono applicati gli stamps.
 * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
 * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo (non usata direttamente).
 * @param prev_solution La soluzione dal passo temporale precedente (non usata direttamente).
 * @param time Il tempo di simulazione corrente (non usata direttamente).
 * @param dt La dimensione del passo temporale (non usata direttamente).
 */
void LDR::getStamps(
    Eigen::MatrixXd& A, Eigen::VectorXd& b,
    const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
    double time, double dt
) {
    // Ottieni gli indici globali per i nodi dell'LDR.
    int idx1 = getNodeIndex(node1);
    int idx2 = getNodeIndex(node2);

    // Controlla la validità degli indici
    if (idx1 == -1 || idx2 == -1) {
        std::cerr << "Errore: Indice di nodo non valido per LDR " << name << std::endl;
        return;
    }

    double current_resistance;

    // Calcola la resistenza attuale dell'LDR in base all'intensità luminosa.
    // Gestisce il caso di luce zero per evitare divisioni per zero o logaritmi di zero.
    if (current_light_intensity <= 1e-9) { // Considera luce quasi zero
        current_resistance = std::numeric_limits<double>::max(); // Resistenza molto alta (circuito aperto)
        std::cout << "LDR " << name << ": Luce quasi zero, resistenza impostata a infinito." << std::endl;
    } else {
        // La formula è R_LDR = R0 * (L / L0)^(-gamma)
        // Equivalente a R_LDR = R0 / (L / L0)^gamma
        current_resistance = R0 * std::pow(current_light_intensity / L0, -gamma);
    }

    // Assicura che la resistenza non sia zero o negativa (può causare problemi numerici)
    if (current_resistance <= 1e-9) { // Limita la resistenza minima
        current_resistance = 1e-9; // 1 nOhm, valore molto piccolo ma non zero
        std::cerr << "Attenzione: LDR " << name << " resistenza calcolata troppo bassa, impostata a " << current_resistance << " Ohm." << std::endl;
    }

    // Calcola la conduttanza
    double G = 1.0 / current_resistance;

    // Applica lo stamp del resistore alla matrice MNA
    // Tra node1 e node2
    A(idx1, idx1) += G;
    A(idx2, idx2) += G;
    A(idx1, idx2) -= G;
    A(idx2, idx1) -= G;
}

/**
 * @brief Imposta l'intensità luminosa attuale dell'LDR.
 *
 * Questo metodo permette di aggiornare l'intensità luminosa esternamente,
 * influenzando la resistenza dell'LDR nel passo di simulazione successivo.
 *
 * @param light_intensity La nuova intensità luminosa in lux. Deve essere >= 0.
 */
void LDR::setLightIntensity(double light_intensity) {
    if (light_intensity < 0.0) {
        std::cerr << "Attenzione: LDR " << name << " tentato di impostare un'intensità luminosa negativa ("
                  << light_intensity << "). Impostazione a 0 lux." << std::endl;
        this->current_light_intensity = 0.0;
    } else {
        this->current_light_intensity = light_intensity;
    }
}
