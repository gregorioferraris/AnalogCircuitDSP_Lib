# circuits/SallenKeyHighPassCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.op_amp import OpAmp

import numpy as np

class SallenKeyHighPassCircuit(Circuit):
    """
    Circuito di un filtro passa-alto attivo di Sallen-Key del secondo ordine.
    Permette di controllare la frequenza di taglio e il fattore Q.
    """
    def __init__(self, name="SallenKey_HP", R1=10000.0, R2=10000.0, C1=1.0e-8, C2=1.0e-8,
                 opamp_gain_resistor_ratio=1.0, opamp_params=None, sample_rate=48000):
        super().__init__(name)
        self.R1_val = float(R1)
        self.R2_val = float(R2)
        self.C1_val = float(C1)
        self.C2_val = float(C2)
        self.opamp_gain_res_ratio = float(opamp_gain_resistor_ratio)

        self.opamp_params = opamp_params if opamp_params is not None else {}
        self.sample_rate = sample_rate

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il filtro Sallen-Key High Pass."""
        self.add_node("Input")
        self.add_node("Output")
        self.add_node("OpAmp_Vout")
        self.add_node("OpAmp_Vplus")
        self.add_node("OpAmp_Vminus")

        # Nodi interni
        self.add_node("C1_R2_Node") # Nodo tra C1 e R2/C2

        # Nodi per il guadagno dell'Op-Amp
        self.add_node("Feedback_Resistor_Node")
        self.add_node("Ground_Resistor_Node")

    def _add_components(self):
        """Aggiunge i componenti (R, C, Op-Amp)."""
        # Capacitor C1
        self.add_component(Capacitor(self.C1_val, sample_rate=self.sample_rate), "Input", "C1_R2_Node")

        # Capacitor C2
        self.add_component(Capacitor(self.C2_val, sample_rate=self.sample_rate), "C1_R2_Node", "OpAmp_Vout")

        # Resistor R1
        self.add_component(Resistor(self.R1_val), "C1_R2_Node", "GND")

        # Resistor R2
        self.add_component(Resistor(self.R2_val), "OpAmp_Vplus", "GND")


        # Op-Amp
        self.op_amp = OpAmp(**self.opamp_params)
        self.add_component(self.op_amp, "OpAmp_Vplus", "OpAmp_Vminus", "OpAmp_Vout")

        # Configurazione del guadagno dell'Op-Amp
        if self.opamp_gain_res_ratio > 1.0 + 1e-9:
            R_g_val = 10000.0
            R_f_val = R_g_val * (self.opamp_gain_res_ratio - 1.0)
            self.add_component(Resistor(R_g_val), "OpAmp_Vminus", "Ground_Resistor_Node")
            self.add_component(Resistor(R_f_val), "Feedback_Resistor_Node", "OpAmp_Vminus")
            self.connect_nodes("Feedback_Resistor_Node", "OpAmp_Vout")
            self.connect_nodes("Ground_Resistor_Node", "GND")
        else: # Guadagno unitario (buffer)
            self.connect_nodes("OpAmp_Vminus", "OpAmp_Vout")


    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        # L'output del filtro è l'uscita dell'Op-Amp
        self.connect_nodes("Output", "OpAmp_Vout")

    def get_input_node(self):
        return self.get_node_id("Input")

    def get_output_node(self):
        return self.get_node_id("Output")

    def calculate_parameters(self):
        """
        Calcola la frequenza di taglio (Fc) e il fattore di merito (Q) del filtro.
        Formule per Sallen-Key passa-alto con guadagno unitario (A=1):
        Fc = 1 / (2 * pi * sqrt(R1*R2*C1*C2))
        Q = sqrt(R1*R2*C1*C2) / (R1 * (C1 + C2))
        """
        if self.R1_val * self.R2_val * self.C1_val * self.C2_val <= 0:
            fc = 0.0
        else:
            fc = 1.0 / (2.0 * np.pi * np.sqrt(self.R1_val * self.R2_val * self.C1_val * self.C2_val))

        A_v = self.opamp_gain_res_ratio
        
        # Per Q, la formula è complessa con guadagno A_v
        # Q = sqrt(R1*R2*C1*C2) / (R1*(C1+C2) + R2*C1*(1-A_v))
        denominator = (self.R1_val * (self.C1_val + self.C2_val)) + (self.R2_val * self.C2_val * (1 - A_v))

        if denominator == 0:
            Q = float('inf')
        else:
            numerator = np.sqrt(self.R1_val * self.R2_val * self.C1_val * self.C2_val)
            Q = numerator / denominator

        return {"cutoff_frequency": fc, "Q": Q, "gain": A_v}

    def set_filter_parameters(self, cutoff_freq, Q_factor, gain=1.0):
        """
        Metodo placeholder per la regolazione dei parametri, come nel passa-basso.
        """
        print(f"DEBUG: Tentativo di impostare Fc={cutoff_freq}Hz, Q={Q_factor}, Gain={gain} (non implementato completamente per ricalcolo R,C)")
        pass
