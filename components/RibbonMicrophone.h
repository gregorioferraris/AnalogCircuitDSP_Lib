// components/RibbonMicrophone.h
#ifndef RIBBON_MICROPHONE_H
#define RIBBON_MICROPHONE_H

#include "Component.h" // Include la classe base Component
#include <functional> // Per std::function per gestire la funzione della pressione sonora
#include <string>

/**
 * @brief Rappresenta un componente Microfono a Nastro in una simulazione di circuito.
 *
 * Questo componente modella l'uscita elettrica di un microfono a nastro come una
 * sorgente di tensione variabile nel tempo (che rappresenta l'input sonoro) in serie con
 * una resistenza di uscita interna e un'induttanza di uscita interna (dovuta principalmente al trasformatore).
 *
 * Ha due terminali elettrici: un nodo di uscita positivo e un nodo di uscita negativo.
 * Per l'analisi transitoria, l'induttore è modellato utilizzando il suo modello di accompagnamento.
 * Questo componente introduce una variabile ausiliaria per la corrente che lo attraversa.
 */
class RibbonMicrophone : public Component {
public:
    /**
     * @brief Costruttore per il componente RibbonMicrophone.
     *
     * Inizializza un microfono a nastro con il suo nome, due nodi di uscita,
     * resistenza di uscita, induttanza di uscita e una funzione che fornisce
     * la pressione sonora istantanea in funzione del tempo.
     *
     * @param name Il nome univoco del microfono.
     * @param positive_node Il nome del nodo di uscita positivo.
     * @param negative_node Il nome del nodo di uscita negativo.
     * @param output_resistance La resistenza di uscita in serie del microfono (Ohm).
     * @param output_inductance L'induttanza di uscita in serie del microfono (Henry).
     * @param sensitivity La sensibilità del microfono (Volt/Pascal).
     * @param sound_pressure_function Una std::function che prende il tempo corrente (double)
     * e restituisce la pressione sonora istantanea (Pascal).
     * @param tolerance_percent Tolleranza opzionale in percentuale (ereditata da Component).
     */
    RibbonMicrophone(const std::string& name,
                     const std::string& positive_node,
                     const std::string& negative_node,
                     double output_resistance,
                     double output_inductance,
                     double sensitivity,
                     std::function<double(double time)> sound_pressure_function,
                     double tolerance_percent = 0.0);

    /**
     * @brief Applica gli "stamps" del RibbonMicrophone alla matrice MNA (A) e al vettore (b).
     *
     * Questo metodo implementa il modello di accompagnamento per la combinazione in serie di una
     * sorgente di tensione variabile nel tempo, una resistenza e un induttore.
     * Aggiunge una variabile ausiliaria per rappresentare la corrente che scorre attraverso il microfono.
     *
     * @param A La matrice MNA a cui vengono applicati gli "stamps".
     * @param b Il vettore lato destro MNA a cui vengono applicati gli "stamps".
     * @param x_current_guess La stima corrente per le tensioni dei nodi e le correnti dei rami.
     * @param prev_solution La soluzione del passo temporale precedente, utilizzata per il modello di accompagnamento dell'induttore.
     * @param time Il tempo di simulazione corrente.
     * @param dt La dimensione del passo temporale.
     */
    void getStamps(
        Eigen::MatrixXd& A, Eigen::VectorXd& b,
        const Eigen::VectorXd& x_current_guess, const Eigen::VectorXd& prev_solution,
        double time, double dt
    ) override;

private:
    double R_out; // Resistenza di uscita (Ohm)
    double L_out; // Induttanza di uscita (Henry)
    double sensitivity; // Sensibilità del microfono (V/Pa)
    // Funzione per fornire la pressione sonora in un dato momento
    std::function<double(double time)> V_sound_pressure_func;

    // Nomi dei nodi per chiarezza (già memorizzati nella classe base Component, ma utili come copie locali)
    std::string pos_node_name;
    std::string neg_node_name;
};

#endif // RIBBON_MICROPHONE_H
