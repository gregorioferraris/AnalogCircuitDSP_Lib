// circuits/InputTransformerCircuit.h
#ifndef INPUT_TRANSFORMER_CIRCUIT_H
#define INPUT_TRANSFORMER_CIRCUIT_H

#include "circuit_solver/subcircuits/Subcircuit.h"
#include "components/Transformer.h" // Usa il componente Transformer effettivo
#include "components/Resistor.h"     // Per qualsiasi resistenza in serie aggiuntiva

/**
 * @class InputTransformerCircuit
 * @brief Modella un circuito di ingresso con un trasformatore.
 *
 * Questo sottocircuito rappresenta un trasformatore ideale con resistenze in serie
 * sugli avvolgimenti primario e secondario. Utilizza il componente `Transformer`
 * per modellare l'accoppiamento magnetico e il rapporto di spire.
 */
class InputTransformerCircuit : public Subcircuit {
public:
    /**
     * @brief Costruttore per il circuito del trasformatore di ingresso.
     * @param name Nome univoco dell'istanza del trasformatore.
     * @param primary_input_node Nome del nodo di ingresso primario esterno.
     * @param primary_ground_node Nome del nodo di massa primario esterno.
     * @param secondary_output_node Nome del nodo di uscita secondario esterno.
     * @param secondary_ground_node Nome del nodo di massa secondario esterno.
     * @param turns_ratio Rapporto di spire (Np/Ns) del trasformatore.
     * @param primary_resistance_ohm Resistenza in serie sul lato primario.
     * @param secondary_resistance_ohm Resistenza in serie sul lato secondario.
     */
    InputTransformerCircuit(const std::string& name,
                            const std::string& primary_input_node,
                            const std::string& primary_ground_node, // Nodo per la massa lato primario
                            const std::string& secondary_output_node,
                            const std::string& secondary_ground_node, // Nodo per la massa lato secondario
                            double turns_ratio,
                            double primary_resistance_ohm, // Resistenza in serie sul lato primario
                            double secondary_resistance_ohm); // Resistenza in serie sul lato secondario

private:
    double turns_ratio_;
    double primary_resistance_ohm_;
    double secondary_resistance_ohm_;

    // Nomi dei nodi interni per i resistori in serie
    std::string primary_res_out_node_;
    std::string secondary_res_out_node_;

    // Puntatore al componente Transformer effettivo
    std::shared_ptr<Transformer> transformer_;

    /**
     * @brief Aggiunge i nodi interni specifici per questo circuito.
     */
    void _add_nodes();

    /**
     * @brief Aggiunge i componenti interni (resistenze e trasformatore).
     */
    void _add_components();

    /**
     * @brief Connette i nodi interni e le porte esterne.
     */
    void _connect_nodes();
};

#endif // INPUT_TRANSFORMER_CIRCUIT_H
