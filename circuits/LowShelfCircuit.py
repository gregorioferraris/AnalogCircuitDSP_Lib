# circuits/LowShelfCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.op_amp import OpAmp

import numpy as np

class LowShelfCircuit(Circuit):
    """
    Circuito di un filtro Low-Shelf attivo.
    Aumenta o diminuisce il guadagno delle frequenze al di sotto di una frequenza di taglio.
    """
    def __init__(self, name="LowShelfFilter", R1=10000.0, R2=10000.0, C1=1.0e-8,
                 opamp_params=None, sample_rate=48000):
        super().__init__(name)
        self.R1_val = float(R1)
        self.R2_val = float(R2)
        self.C1_val = float(C1)
        self.opamp_params = opamp_params if opamp_params is not None else {}
        self.sample_rate = sample_rate

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il filtro Low-Shelf."""
        self.add_node("Input")
        self.add_node("Output")
        self.add_node("OpAmp_Vout")
        self.add_node("OpAmp_Vminus")
        self.add_node("OpAmp_Vplus")

        # Nodi interni
        self.add_node("Feedback_Node_1") # Tra R1 e C1
        self.add_node("Feedback_Node_2") # Tra C1 e R2

    def _add_components(self):
        """Aggiunge i componenti (R, C, Op-Amp)."""
        # Op-Amp
        self.op_amp = OpAmp(**self.opamp_params)
        self.add_component(self.op_amp, "OpAmp_Vplus", "OpAmp_Vminus", "OpAmp_Vout")

        # Assumiamo una configurazione di Op-Amp invertente per il filtro shelf.
        # R1 in ingresso
        self.add_component(Resistor(self.R1_val), "Input", "OpAmp_Vminus")

        # C1 e R2 in parallelo in feedback
        self.add_component(Capacitor(self.C1_val, sample_rate=self.sample_rate), "OpAmp_Vout", "Feedback_Node_1")
        self.add_component(Resistor(self.R2_val), "OpAmp_Vout", "Feedback_Node_1")
        self.connect_nodes("Feedback_Node_1", "OpAmp_Vminus") # Il parallelo va a Vminus

        # L'ingresso non invertente a GND
        self.connect_nodes("OpAmp_Vplus", "GND")

    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        self.connect_nodes("Output", "OpAmp_Vout")

    def get_input_node(self):
        return self.get_node_id("Input")

    def get_output_node(self):
        return self.get_node_id("Output")

    def calculate_parameters(self):
        """
        Calcola la frequenza di taglio (Fc) e il guadagno del Low-Shelf.
        Per un Low-Shelf con feedback parallelo C1, R2 e ingresso R1:
        Gain_low_freq = -R2 / R1
        Gain_high_freq = - (R2 || C1) / R1 (R2 parallelo C1) -> -R2 / R1 (quando C1 è un corto)
        Fc = 1 / (2 * pi * R2 * C1)
        """
        # Guadagno a basse frequenze (quando C1 è un aperto):
        if self.R1_val == 0:
            gain_low_freq_linear = float('inf')
        else:
            gain_low_freq_linear = -self.R2_val / self.R1_val
        gain_low_freq_db = 20 * np.log10(abs(gain_low_freq_linear))

        # Frequenza di taglio (Fc dove la pendenza inizia/finisce)
        if self.R2_val * self.C1_val == 0:
            fc = 0.0
        else:
            fc = 1.0 / (2.0 * np.pi * self.R2_val * self.C1_val)

        return {"cutoff_frequency": fc, "gain_db": gain_low_freq_db}

    def set_parameters(self, cutoff_freq, gain_db=0.0):
        """
        Metodo per sintetizzare R e C per ottenere Fc e Gain desiderati.
        """
        print(f"DEBUG: Sintesi di LowShelf per Fc={cutoff_freq}Hz, Gain={gain_db}dB (non implementato)")
        # Esempio di sintesi:
        # C_ref = 1.0e-8 # Fissa C di riferimento
        # R2_val = 1.0 / (2.0 * np.pi * cutoff_freq * C_ref)
        # R1_val = R2_val / abs(10**(gain_db/20.0))
        # Aggiornare self.R1_val, self.R2_val, self.C1_val.
        pass
