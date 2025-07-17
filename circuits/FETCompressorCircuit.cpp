// circuits/FETCompressorCircuit.cpp
#include "FETCompressorCircuit.h"
#include <iostream>

FETCompressorCircuit::FETCompressorCircuit(
    const std::string& name,
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
    double sample_rate)
    : Subcircuit(name, audio_input_node, audio_output_node, ground_node), // Le porte esterne sono input, output, ground
      mosfet_vt_(mosfet_vt), mosfet_kn_(mosfet_kn), mosfet_lambda_val_(mosfet_lambda_val),
      sidechain_diode_is_(sidechain_diode_is), sidechain_diode_n_(sidechain_diode_n),
      sidechain_attack_ms_(sidechain_attack_ms), sidechain_release_ms_(sidechain_release_ms),
      sidechain_filter_cap_ref_(sidechain_filter_cap_ref),
      vca_input_resistor_(vca_input_resistor),
      vca_feedback_resistor_fixed_(vca_feedback_resistor_fixed),
      gate_bias_resistor_gnd_(gate_bias_resistor_gnd),
      gate_bias_resistor_vcc_(gate_bias_resistor_vcc),
      bias_voltage_vcc_(bias_voltage_vcc),
      sample_rate_(sample_rate)
{
    std::cout << "Inizializzazione FETCompressorCircuit: " << name << std::endl;
    _add_nodes();
    _add_components();
    _connect_nodes();
}

void FETCompressorCircuit::_add_nodes() {
    // Nodi per la side-chain
    // sc_audio_input_node_ è concettualmente l'input della sidechain, che sarà collegato all'input audio principale.
    // Non è un nodo *aggiuntivo* ma un alias per la porta di input del sottocircuito sidechain.
    sc_control_voltage_output_node_ = getName() + "_SC_Control_Voltage_Output";

    // Nodi per il MOSFET (Drain, Gate, Source)
    mosfet_drain_node_ = getName() + "_MOSFET_Drain";
    mosfet_gate_node_ = getName() + "_MOSFET_Gate";
    mosfet_source_node_ = getName() + "_MOSFET_Source";

    // Nodi per il circuito VCA (Op-Amp e resistori)
    vca_opamp_vminus_node_ = getName() + "_VCA_OpAmp_Vminus";
    vca_opamp_vplus_node_ = getName() + "_VCA_OpAmp_Vplus";
    vca_opamp_vout_node_ = getName() + "_VCA_OpAmp_Vout";

    // Nodo interno per la rete di bias del gate
    gate_bias_mid_node_ = getName() + "_Gate_Bias_Mid";
}

void FETCompressorCircuit::_add_components() {
    // 1. Side-chain (Diode Bridge Compressor Sidechain)
    // L'input della sidechain è il nodo di input audio esterno del compressore (getPortNames()[0]).
    sidechain_circuit_ = std::make_shared<DiodeBridgeCompressorSidechain>(
        getName() + "_Sidechain",
        getPortNames()[0], // Ingresso audio principale del compressore
        sc_control_voltage_output_node_, // Nodo interno che sarà l'output di controllo della side-chain
        getPortNames()[2], // Massa esterna
        sidechain_diode_is_, sidechain_diode_n_,
        sidechain_attack_ms_, sidechain_release_ms_,
        sidechain_filter_cap_ref_, sample_rate_);
    // Aggiungi il sottocircuito della side-chain al circuito interno di questo compressore
    internalCircuit.addSubcircuit(*sidechain_circuit_);

    // 2. MOSFET (che agisce come resistore variabile)
    // Tipo "NMOS" per simulare il comportamento di un JFET a canale N con Vt negativo.
    // Il nodo Bulk è collegato al Source per semplicità (comune per un comportamento simile a JFET).
    mosfet_ = std::make_shared<MOSFET>(
        getName() + "_MOSFET",
        mosfet_drain_node_, mosfet_gate_node_, mosfet_source_node_, mosfet_source_node_, // Drain, Gate, Source, Bulk (collegato al Source)
        "NMOS", mosfet_kn_, mosfet_vt_, mosfet_lambda_val_,
        0.0, 0.6, 100e-6, 100e-6); // Gamma, Phi, W, L, Tolerance di default
    addInternalComponent(mosfet_);

    // 3. Op-Amp VCA e resistori associati
    vca_opamp_ = std::make_shared<OpAmp>(
        getName() + "_VCA_OpAmp",
        vca_opamp_vplus_node_, vca_opamp_vminus_node_, vca_opamp_vout_node_);
    addInternalComponent(vca_opamp_);

    // Resistore di ingresso per l'Op-Amp VCA
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_VCA_R_in", getPortNames()[0], vca_opamp_vminus_node_, vca_input_resistor_));

    // Resistore di feedback fisso in serie con il MOSFET
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_VCA_R_fb_fixed", vca_opamp_vout_node_, mosfet_drain_node_, vca_feedback_resistor_fixed_));

    // Rete di bias del gate per il MOSFET
    // Resistore da Gate_Bias_Mid a Massa
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_Gate_Bias_R_GND", gate_bias_mid_node_, getPortNames()[2], gate_bias_resistor_gnd_));
    // Resistore da Gate_Bias_Mid a Tensione di Bias (VCC)
    addInternalComponent(std::make_shared<Resistor>(
        getName() + "_Gate_Bias_R_VCC", gate_bias_mid_node_, getName() + "_Bias_VCC_Node", gate_bias_resistor_vcc_));

    // Sorgente di tensione di bias (DC)
    bias_source_ = std::make_shared<VoltageSource>(
        getName() + "_Bias_VCC_Source", getName() + "_Bias_VCC_Node", getPortNames()[2],
        VoltageSource::DC, bias_voltage_vcc_, 0.0, 0.0, sample_rate_);
    addInternalComponent(bias_source_);
}

void FETCompressorCircuit::_connect_nodes() {
    // Connetti l'ingresso non invertente dell'Op-Amp VCA a massa
    internalCircuit.connectNodes(vca_opamp_vplus_node_, getPortNames()[2]);

    // Connetti il Source del MOSFET all'ingresso invertente dell'Op-Amp VCA
    internalCircuit.connectNodes(mosfet_source_node_, vca_opamp_vminus_node_);

    // Connetti l'uscita della tensione di controllo della side-chain alla rete di bias del gate del MOSFET
    internalCircuit.connectNodes(sc_control_voltage_output_node_, mosfet_gate_node_); // Collegato direttamente al gate

    // Mappa il nodo di uscita interno dell'Op-Amp VCA alla porta di uscita audio esterna del compressore
    internalCircuit.connectNodes(vca_opamp_vout_node_, getPortNames()[1]);
}
