# circuits/FETCompressorCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.op_amp import OpAmp # Se usiamo Op-Amp nel percorso audio
from components.mosfet import MOSFET # Il TUO componente MOSFET aggiornato
from circuits.DiodeBridgeCompressorSidechain import DiodeBridgeCompressorSidechain

import numpy as np

class FETCompressorCircuit(Circuit):
    """
    Circuito di un compressore audio basato su FET.
    Incorpora un MOSFET (configurabile come JFET a impoverimento) come elemento di guadagno variabile
    e una side-chain per generare la tensione di controllo.

    Componenti chiave:
    - MOSFET: agisce come resistore variabile (VCR) nel percorso del segnale audio.
    - DiodeBridgeCompressorSidechain: rileva l'inviluppo del segnale per la tensione di controllo.
    - Rete di controllo: collega la tensione di inviluppo al gate del MOSFET.
    - Circuito VCA: in questo caso, il MOSFET stesso inserito in una configurazione di VCA.
    """
    def __init__(self, name="FET_Compressor",
                 # Parametri del MOSFET (regolali per comportamento JFET)
                 mosfet_vt=-3.0, # Tipico per JFET N-channel a impoverimento (es. J201, 2N5457)
                 mosfet_kn=0.005, # Transconduttanza
                 mosfet_lambda_val=0.01, # Modulazione di canale, spesso piccola
                 
                 # Parametri della side-chain (ponte diodi + RC)
                 sidechain_diode_is=1e-14, sidechain_diode_n=1.0,
                 sidechain_attack_ms=5.0, # Tempi rapidi per FET
                 sidechain_release_ms=50.0, # Tempi rapidi per FET
                 sidechain_filter_cap_ref=1.0e-6, # Capacità di riferimento per il filtro RC della side-chain

                 # Parametri del circuito VCA (con il MOSFET)
                 vca_input_resistor=10000.0, # Resistenza di ingresso per l'Op-Amp invertente
                 vca_feedback_resistor_fixed=20000.0, # Resistenza fissa in serie al MOSFET nel feedback
                 # Aggiungiamo un resistore per l'offset DC sul gate del FET, per il "Threshold"
                 gate_bias_resistor_gnd=100000.0, # Resistenza dal gate a GND
                 gate_bias_resistor_vcc=100000.0, # Resistenza dal gate a VCC (o sorgente di bias)
                 bias_voltage_vcc=5.0, # Tensione di bias (per polarizzare il FET)

                 sample_rate=48000):
        super().__init__(name)
        self.sample_rate = sample_rate

        self.mosfet_params = {
            "Vt": mosfet_vt,
            "Kn": mosfet_kn,
            "lambda_val": mosfet_lambda_val
        }
        
        self.sidechain_params = {
            "diode_is": sidechain_diode_is,
            "diode_n": sidechain_diode_n,
            "sample_rate": sample_rate,
            "filter_capacitor_value": sidechain_filter_cap_ref # Capacità per la side-chain
        }
        self.sidechain_attack_ms = sidechain_attack_ms
        self.sidechain_release_ms = sidechain_release_ms

        self.vca_R_in_val = float(vca_input_resistor)
        self.vca_R_fb_fixed_val = float(vca_feedback_resistor_fixed)
        self.gate_bias_R_gnd_val = float(gate_bias_resistor_gnd)
        self.gate_bias_R_vcc_val = float(gate_bias_resistor_vcc)
        self.bias_voltage_vcc = float(bias_voltage_vcc)

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()
        
        # Setup iniziale dei tempi della side-chain
        # Questo creerà l'istanza della sidechain con i R/C calcolati
        self.set_attack_release_times(self.sidechain_attack_ms, self.sidechain_release_ms)

    def _add_nodes(self):
        """Aggiunge i nodi per l'intero compressore."""
        self.add_node("Audio_Input")
        self.add_node("Audio_Output")

        # Nodi per la side-chain
        self.add_node("SC_Audio_Input") # Ingresso del segnale audio per la side-chain
        self.add_node("SC_Control_Voltage_Output") # Uscita della tensione di controllo DC dalla side-chain

        # Nodi per il MOSFET (Drain, Gate, Source)
        self.add_node("MOSFET_Drain")
        self.add_node("MOSFET_Gate")
        self.add_node("MOSFET_Source")

        # Nodi per il circuito VCA (Op-Amp e resistori)
        self.add_node("VCA_OpAmp_Vminus") # Ingresso invertente dell'Op-Amp VCA
        self.add_node("VCA_OpAmp_Vplus")  #
