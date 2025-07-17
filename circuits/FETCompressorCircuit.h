// circuits/FETCompressorCircuit.h
#ifndef FET_COMPRESSOR_CIRCUIT_H
#define FET_COMPRESSOR_CIRCUIT_H

#include "circuit_solver/subcircuits/Subcircuit.h"
#include "components/Resistor.h"
#include "components/Capacitor.h"
#include "components/OpAmp.h"
#include "components/MOSFET.h"
#include "components/VoltageSource.h" // Per la tensione di bias
#include "DiodeBridgeCompressorSidechain.h" // Includi il sottocircuito della side-chain

/**
 * @class FETCompressorCircuit
 * @brief Modella un circuito di compressore audio basato su FET.
 *
 * Questo circuito incorpora un MOSFET come elemento a guadagno variabile (VCA) e una side-chain
 * per generare la tensione di controllo. La side-chain rileva l'inviluppo del segnale audio
 * e la sua uscita controlla la conduttanza del MOSFET, influenzando così il guadagno del VCA.
 */
class FETCompressorCircuit : public Subcircuit {
public:
    /**
     * @brief Costruttore per il circuito del compressore FET.
     * @param name Nome univoco dell'istanza del compressore.
     * @param audio_input_node Nome del nodo di ingresso audio esterno.
     * @param audio_output_node Nome del nodo di uscita audio esterno.
     * @param ground_node Nome del nodo di massa esterno.
     * @param mosfet_vt Tensione di soglia del MOSFET (Vt).
     * @param mosfet_kn Parametro di transconduttanza del MOSFET (Kn).
     * @param mosfet_lambda_val Parametro di modulazione della lunghezza del canale (Lambda).
     * @param sidechain_diode_is Corrente di saturazione dei diodi della side-chain.
     * @param sidechain_diode_n Fattore di idealità dei diodi della side-chain.
     * @param sidechain_attack_ms Tempo di attacco della side-chain in ms.
     * @param sidechain_release_ms Tempo di rilascio della side-chain in ms.
     * @param sidechain_filter_cap_ref Capacità di riferimento per il filtro della side-chain.
     * @param vca_input_resistor Valore della resistenza di ingresso del VCA.
     * @param vca_feedback_resistor_fixed Valore della resistenza fissa nel feedback del VCA.
     * @param gate_bias_resistor_gnd Valore della resistenza di bias del gate verso massa.
     * @param gate_bias_resistor_vcc Valore della resistenza di bias del gate verso VCC.
     * @param bias_voltage_vcc Tensione di bias VCC.
     * @param sample_rate Frequenza di campionamento audio in Hz.
     */
    FETCompressorCircuit(const std::string& name,
                         const std::string& audio_input_node,
                         const std::string& audio_output_node,
                         const std::string& ground_node,
                         double mosfet_vt, double mosfet_kn, double mosfet_lambda_val,
                         double sidechain_diode_is, double sidechain_diode_n,
                         double sidechain_attack_ms, double sidechain_release_ms,
                         double sidechain_filter_cap_ref,
                         double vca_input_resistor,
                         double vca_feedback_resistor_fixed,
                         double gate_bias_resistor_gnd,
                         double gate_bias_resistor_vcc,
                         double bias_voltage_vcc,
                         double sample_rate);

private:
    // Parametri del circuito
    double mosfet_vt_;
    double mosfet_kn_;
    double mosfet_lambda_val_;
    double sidechain_diode_is_;
    double sidechain_diode_n_;
    double sidechain_attack_ms_;
    double sidechain_release_ms_;
    double sidechain_filter_cap_ref_;
    double vca_input_resistor_;
    double vca_feedback_resistor_fixed_;
    double gate_bias_resistor_gnd_;
    double gate_bias_resistor_vcc_;
    double bias_voltage_vcc_;
    double sample_rate_;

    // Nomi dei nodi interni
    std::string sc_audio_input_node_; // Ingresso audio per la side-chain (mappato all'input principale)
    std::string sc_control_voltage_output_node_; // Uscita della tensione di controllo dalla side-chain
    std::string mosfet_drain_node_;
    std::string mosfet_gate_node_;
    std::string mosfet_source_node_;
    std::string vca_opamp_vminus_node_;
    std::string vca_opamp_vplus_node_;
    std::string vca_opamp_vout_node_;
    std::string gate_bias_mid_node_; // Nodo intermedio per la rete di bias del gate

    // Puntatori ai sottocircuiti e componenti interni (opzionali, ma utili per accesso futuro)
    std::shared_ptr<DiodeBridgeCompressorSidechain> sidechain_circuit_;
    std::shared_ptr<MOSFET> mosfet_;
    std::shared_ptr<OpAmp> vca_opamp_;
    std::shared_ptr<VoltageSource> bias_source_;

    /**
     * @brief Aggiunge i nodi interni specifici per questo circuito.
     */
    void _add_nodes();

    /**
     * @brief Aggiunge i componenti e i sottocircuiti interni.
     */
    void _add_components();

    /**
     * @brief Connette i nodi interni e le porte esterne.
     */
    void _connect_nodes();
};

#endif // FET_COMPRESSOR_CIRCUIT_H
