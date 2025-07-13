# circuits/TwinTNotchCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.op_amp import OpAmp

import numpy as np

class TwinTNotchCircuit(Circuit):
    """
    Circuito di un filtro Notch (elimina-banda) basato su una rete Twin-T passiva
    in combinazione con un amplificatore operazionale.
    La rete Twin-T è una configurazione di R e C che crea una null alla frequenza centrale.
    L'Op-Amp può migliorare il Q (selettività) e la profondità del notch.
    """
    def __init__(self, name="TwinT_Notch", R=10000.0, C=1.0e-8, Q_boost_resistor_ratio=0.0,
                 opamp_params=None, sample_rate=48000):
        super().__init__(name)
        self.R_val = float(R)
        self.C_val = float(C)
        self.Q_boost_ratio = float(Q_boost_resistor_ratio) # Rapporto per aumentare il Q
        self.opamp_params = opamp_params if opamp_params is not None else {}
        self.sample_rate = sample_rate

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il filtro Twin-T."""
        self.add_node("Input")
        self.add_node("Output")
        self.add_node("OpAmp_Vout")
        self.add_node("OpAmp_Vminus") # Ingresso invertente (dove si applica il Twin-T)
        self.add_node("OpAmp_Vplus")  # Ingresso non invertente (a GND)

        # Nodi interni della rete Twin-T
        self.add_node("T1_Node1")
        self.add_node("T1_Node2")
        self.add_node("T2_Node1")

        # Nodi per il Q-boost (se usato)
        self.add_node("Q_Boost_Resistor_Node")

    def _add_components(self):
        """Aggiunge i componenti della rete Twin-T e l'Op-Amp."""
        # --- Parte della prima "T" (Resistori) ---
        self.add_component(Resistor(self.R_val), "Input", "T1_Node1")
        self.add_component(Resistor(self.R_val), "T1_Node1", "T1_Node2")
        self.add_component(Resistor(self.R_val / 2.0), "T1_Node2", "GND") # Questo è R/2

        # --- Parte della seconda "T" (Condensatori) ---
        self.add_component(Capacitor(self.C_val, sample_rate=self.sample_rate), "Input", "T2_Node1")
        self.add_component(Capacitor(self.C_val, sample_rate=self.sample_rate), "T2_Node1", "OpAmp_Vminus") # L'uscita del Twin-T va all'OpAmp Vminus
        self.add_component(Capacitor(self.C_val * 2.0, sample_rate=self.sample_rate), "T1_Node1", "T2_Node1") # Questo è 2C, tra i nodi comuni

        # Op-Amp
        self.op_amp = OpAmp(**self.opamp_params)
        self.add_component(self.op_amp, "OpAmp_Vplus", "OpAmp_Vminus", "OpAmp_Vout")

        # Feedback dell'Op-Amp per il guadagno e il Q-boost
        # Se Q_boost_ratio = 0, l'Op-Amp è un follower di tensione.
        # Se Q_boost_ratio > 0, si applica un feedback positivo (o negativo a seconda della configurazione)
        # per aumentare il Q.
        # Per un Twin-T attivo, il feedback è solitamente dall'uscita dell'Op-Amp a T1_Node2 (la giunzione R/2 con GND).
        if self.Q_boost_ratio > 1e-9: # Se vogliamo aumentare il Q
            # Riferimento: R_q_boost = R / (4 * Q_boost_ratio)
            # Questo è per un Q-boost con feedback.
            # Una configurazione comune è un resistore tra OpAmp_Vout e T1_Node2.
            # O un partitore resistivo per creare il Q-boost.
            # Aggiungiamo un resistore di feedback per aumentare il Q (guadagno del filtro).
            R_feedback_val = self.R_val * (self.Q_boost_ratio) # Semplificazione, necessita calcolo preciso
            self.add_component(Resistor(R_feedback_val), "OpAmp_Vout", "T1_Node2")


    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        self.connect_nodes("OpAmp_Vplus", "GND") # Ingresso non invertente a GND
        self.connect_nodes("Output", "OpAmp_Vout")

    def calculate_parameters(self):
        """
        Calcola la frequenza centrale (Fc) e il fattore di merito (Q) del notch.
        Per la rete Twin-T, Fc = 1 / (2 * pi * R * C) quando R_top=R_mid=R, R_bottom=R/2, C_top=C_mid=C, C_bottom=2C.
        """
        if self.R_val * self.C_val <= 0:
            fc = 0.0
        else:
            fc = 1.0 / (2.0 * np.pi * self.R_val * self.C_val)

        # Il Q in un Twin-T attivo dipende dalla quantità di feedback.
        # Per un Q-boost, il rapporto di Q_boost_ratio influenza il Q.
        # Q = f(Q_boost_ratio)
        return {"center_frequency": fc, "Q": "Determined by R, C, and Q_boost_ratio", "depth_db": "Determined by R, C, and Q_boost_ratio"}

    def set_parameters(self, center_freq, Q_factor):
        """
        Metodo per sintetizzare R e C e il Q_boost_ratio per ottenere Fc e Q.
        """
        print(f"DEBUG: Sintesi di TwinT Notch per Fc={center_freq}Hz, Q={Q_factor} (non implementato)")
        # Esempio di sintesi:
        # C_ref = 1.0e-8 # Fissa C di riferimento
        # R_val = 1.0 / (2.0 * np.pi * center_freq * C_ref)
        # Aggiornare self.R_val, self.C_val.
        # E il Q_boost_ratio per il Q desiderato.
        pass
