# circuits/SallenKeyBandPassCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.op_amp import OpAmp

import numpy as np

class SallenKeyBandPassCircuit(Circuit):
    """
    Circuito di un filtro passa-banda attivo di tipo Sallen-Key.
    """
    def __init__(self, name="SallenKey_BandPass", R1=10000.0, R2=10000.0, R3=10000.0,
                 C1=1.0e-8, C2=1.0e-8, opamp_params=None, sample_rate=48000):
        super().__init__(name)
        self.R1_val = float(R1)
        self.R2_val = float(R2)
        self.R3_val = float(R3)
        self.C1_val = float(C1)
        self.C2_val = float(C2)
        self.opamp_params = opamp_params if opamp_params is not None else {}
        self.sample_rate = sample_rate

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il filtro Sallen-Key Band Pass."""
        self.add_node("Input")
        self.add_node("Output")
        self.add_node("OpAmp_Vout")
        self.add_node("OpAmp_Vplus")
        self.add_node("OpAmp_Vminus")

        # Nodi interni
        self.add_node("Node_A") # Tra C1, R1
        self.add_node("Node_B") # Tra C2, R2

    def _add_components(self):
        """Aggiunge i componenti (R, C, Op-Amp)."""
        # C1 in ingresso
        self.add_component(Capacitor(self.C1_val, sample_rate=self.sample_rate), "Input", "Node_A")

        # R1 in serie
        self.add_component(Resistor(self.R1_val), "Node_A", "OpAmp_Vplus")

        # C2 a GND
        self.add_component(Capacitor(self.C2_val, sample_rate=self.sample_rate), "Node_A", "GND")

        # R2 in feedback da Vout a Vminus
        self.add_component(Resistor(self.R2_val), "OpAmp_Vout", "OpAmp_Vminus")

        # R3 da Vminus a GND (o riferimento)
        self.add_component(Resistor(self.R3_val), "OpAmp_Vminus", "GND")


        # Op-Amp
        self.op_amp = OpAmp(**self.opamp_params)
        self.add_component(self.op_amp, "OpAmp_Vplus", "OpAmp_Vminus", "OpAmp_Vout")

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
        Calcola la frequenza centrale (Fc), il fattore di merito (Q) e il guadagno.
        Formule per Sallen-Key Band-Pass. Anche qui, le formule sono complesse
        e richiedono una sintesi per ottenere Fc, Q, Gain da R, C.
        """
        # Semplificazione: assumendo R1=R3=R e C1=C2=C
        # Fc = 1 / (2 * pi * R * C)
        # Q = 0.5
        # Gain = 1
        
        # Una formula generica per Fc (assumendo il guadagno unitario dell'OpAmp per semplicità):
        fc = 1.0 / (2.0 * np.pi * np.sqrt(self.R1_val * self.R2_val * self.C1_val * self.C2_val))
        
        # Il Q e il guadagno dipendono fortemente dalle relazioni tra R1,R2,R3,C1,C2
        # e dal guadagno dell'OpAmp. Questo modello specifico di Sallen-Key Passa-Banda
        # ha formule complesse.
        # Useremo un placeholder per Q e Gain.
        return {"center_frequency": fc, "Q": "Determined by R, C values", "gain_db": "Determined by R, C values"}

    def set_parameters(self, center_freq, Q_factor, gain_db=0.0):
        """
        Metodo placeholder per la regolazione dei parametri, come negli altri filtri.
        """
        print(f"DEBUG: Sintesi di SallenKey BandPass per Fc={center_freq}Hz, Q={Q_factor}, Gain={gain_db}dB (non implementato)")
        pass
