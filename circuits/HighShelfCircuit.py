# circuits/HighShelfCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.op_amp import OpAmp

import numpy as np

class HighShelfCircuit(Circuit):
    """
    Circuito di un filtro High-Shelf attivo.
    Aumenta o diminuisce il guadagno delle frequenze al di sopra di una frequenza di taglio.
    """
    def __init__(self, name="HighShelfFilter", R1=10000.0, R2=10000.0, C1=1.0e-8,
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
        """Aggiunge i nodi specifici per il filtro High-Shelf."""
        self.add_node("Input")
        self.add_node("Output")
        self.add_node("OpAmp_Vout")
        self.add_node("OpAmp_Vminus")
        self.add_node("OpAmp_Vplus")

        # Nodi interni
        self.add_node("Feedback_Node") # Tra R1, C1, R2

    def _add_components(self):
        """Aggiunge i componenti (R, C, Op-Amp)."""
        # Op-Amp
        self.op_amp = OpAmp(**self.opamp_params)
        self.add_component(self.op_amp, "OpAmp_Vplus", "OpAmp_Vminus", "OpAmp_Vout")

        # Assumiamo una configurazione di Op-Amp invertente per il filtro shelf.
        # Ingresso con R1 e C1 in serie
        self.add_component(Resistor(self.R1_val), "Input", "Feedback_Node")
        self.add_component(Capacitor(self.C1_val, sample_rate=self.sample_rate), "Feedback_Node", "OpAmp_Vminus")

        # R2 in feedback da Vout a Vminus
        self.add_component(Resistor(self.R2_val), "OpAmp_Vout", "OpAmp_Vminus")

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
        Calcola la frequenza di taglio (Fc) e il guadagno del High-Shelf.
        Per un High-Shelf con feedback R2 e ingresso R1 in serie con C1:
        Gain_high_freq = -R2 / R1
        Gain_low_freq = -R2 / (R1 + Z_C1) -> -R2/R1 (se C1 è un aperto)
        Fc = 1 / (2 * pi * R1 * C1)
        """
        # Guadagno ad alte frequenze (quando C1 è un corto):
        if self.R1_val == 0:
            gain_high_freq_linear = float('inf')
        else:
            gain_high_freq_linear = -self.R2_val / self.R1_val
        gain_high_freq_db = 20 * np.log10(abs(gain_high_freq_linear))

        # Frequenza di taglio (Fc dove la pendenza inizia/finisce)
        if self.R1_val * self.C1_val == 0:
            fc = 0.0
        else:
            fc = 1.0 / (2.0 * np.pi * self.R1_val * self.C1_val)

        return {"cutoff_frequency": fc, "gain_db": gain_high_freq_db}

    def set_parameters(self, cutoff_freq, gain_db=0.0):
        """
        Metodo per sintetizzare R e C per ottenere Fc e Gain desiderati.
        """
        print(f"DEBUG: Sintesi di HighShelf per Fc={cutoff_freq}Hz, Gain={gain_db}dB (non implementato)")
        # Esempio di sintesi:
        # C_ref = 1.0e-8 # Fissa C di riferimento
        # R1_val = 1.0 / (2.0 * np.pi * cutoff_freq * C_ref)
        # R2_val = R1_val * abs(10**(gain_db/20.0))
        # Aggiornare self.R1_val, self.R2_val, self.C1_val.
        pass
