// components/DynamicMicrophone.h
#ifndef DYNAMIC_MICROPHONE_H
#define DYNAMIC_MICROPHONE_H

#include "Component.h"   // Include la classe base Component
#include <Eigen/Dense>   // Per le operazioni con matrici e vettori Eigen
#include <string>        // Per std::string
#include <functional>    // Per std::function

/**
 * @class DynamicMicrophone
 * @brief Rappresenta un microfono dinamico come una sorgente di tensione variabile nel tempo.
 *
 * Questo componente modella un microfono dinamico come una sorgente di tensione ideale
 * la cui tensione di uscita è una funzione del tempo, simulando la conversione
 * della pressione sonora in un segnale elettrico. Introduce una variabile ausiliaria
 * per la sua corrente di ramo nella matrice MNA.
 */
class DynamicMicrophone : public Component {
public:
    // Funzione per definire la tensione di uscita in funzione del tempo.
    // Questo simula l'input sonoro.
    std::function<double(double)> voltage_function;

    /**
     * @brief Costruttore per il componente DynamicMicrophone.
     * @param name Il nome univoco del microfono.
     * @param node1 Il nome del nodo del terminale positivo.
     * @param node2 Il nome del nodo del terminale negativo.
     * @param v_func Una funzione che restituisce la tensione di uscita del microfono a un dato tempo.
     * @param tolerance_percent Percentuale di tolleranza opzionale (non usata direttamente per la sorgente ideale).
     */
    DynamicMicrophone(const std::string& name, const std::string& node1, const std::string& node2,
                      std::function<double(double)> v_func, double tolerance_percent = 0.0);

    /**
     * @brief Crea una copia profonda dell'oggetto DynamicMicrophone.
     * @return Un puntatore a un nuovo oggetto DynamicMicrophone, che è una copia dell'istanza corrente.
     */
    Component* clone() const override { return new DynamicMicrophone(*this); }

    /**
     * @brief Applica gli "stamps" del microfono dinamico (come sorgente di tensione)
     * alla matrice MNA (A) e al vettore (b).
     *
     * Questo componente introduce una variabile ausiliaria per la sua corrente di ramo.
     * Gli stamps impongono il vincolo di tensione $V(node1) - V(node2) = \text{voltage\_function(time)}$.
     *
     * @param A La matrice MNA a cui vengono applicati gli stamps.
     * @param b Il vettore lato destro MNA a cui vengono applicati gli stamps.
     * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti di ramo (non usata direttamente qui).
     * @param prev_solution La soluzione dal passo temporale precedente (non usata qui).
     * @param time Il tempo di simulazione corrente.
     * @param dt La dimensione del passo temporale (non usata qui).
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

    /**
     * @brief Aggiorna lo stato interno del microfono dinamico.
     *
     * Essendo una sorgente di tensione ideale, l'uscita del microfono è determinata unicamente
     * dalla funzione di tensione fornita in funzione del tempo e non ha uno stato interno
     * che evolve in base alla propria tensione o corrente.
     * Questo metodo è lasciato vuoto.
     *
     * @param v_curr La tensione corrente attraverso il componente (non usata).
     * @param i_curr La corrente corrente che attraversa il componente (non usata).
     */
    void updateState(double v_curr, double i_curr) override;
};

#endif // DYNAMIC_MICROPHONE_H
