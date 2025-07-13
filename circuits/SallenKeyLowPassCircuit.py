# circuits/SallenKeyLowPassCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.op_amp import OpAmp

import numpy as np

class SallenKeyLowPassCircuit(Circuit):
    """
    Circuito di un filtro passa-basso attivo di Sallen-Key del secondo ordine.
    Permette di controllare la frequenza di taglio e il fattore Q (damping/risonanza).
    """
    def __init__(self, name="SallenKey_LP", R1=10000.0, R2=10000.0, C1=1.0e-8, C2=1.0e-8,
                 opamp_gain_resistor_ratio=1.0, opamp_params=None, sample_rate=48000):
        super().__init__(name)
        self.R1_val = float(R1)
        self.R2_val = float(R2)
        self.C1_val = float(C1)
        self.C2_val = float(C2)
        self.opamp_gain_res_ratio = float(opamp_gain_resistor_ratio) # Per guadagno non unitario

        self.opamp_params = opamp_params if opamp_params is not None else {}
        self.sample_rate = sample_rate

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il filtro Sallen-Key."""
        self.add_node("Input")
        self.add_node("Output")
        self.add_node("OpAmp_Vout") # Uscita dell'Op-Amp
        self.add_node("OpAmp_Vplus") # Ingresso non invertente
        self.add_node("OpAmp_Vminus") # Ingresso invertente

        # Nodi interni
        self.add_node("R1_C1_Node") # Nodo tra R1 e C1/R2

        # Nodi per il guadagno dell'Op-Amp (se non unitario)
        self.add_node("Feedback_Resistor_Node")
        self.add_node("Ground_Resistor_Node")

        # GND e VCC/VEE (per l'Op-Amp, se necessario un modello più complesso) sono impliciti
        # Qui l'Op-Amp ideale non ha bisogno di alimentazione esplicita per i nodi, solo Vsat.

    def _add_components(self):
        """Aggiunge i componenti (R, C, Op-Amp)."""
        # Resistor R1
        self.add_component(Resistor(self.R1_val), "Input", "R1_C1_Node")

        # Resistor R2
        self.add_component(Resistor(self.R2_val), "R1_C1_Node", "OpAmp_Vplus")

        # Capacitor C1
        self.add_component(Capacitor(self.C1_val, sample_rate=self.sample_rate), "R1_C1_Node", "OpAmp_Vout")

        # Capacitor C2
        self.add_component(Capacitor(self.C2_val, sample_rate=self.sample_rate), "OpAmp_Vplus", "GND")

        # Op-Amp
        self.op_amp = OpAmp(**self.opamp_params)
        self.add_component(self.op_amp, "OpAmp_Vplus", "OpAmp_Vminus", "OpAmp_Vout")

        # Implementazione di un buffer o di un guadagno non unitario con l'Op-Amp
        # Per un Sallen-Key standard, l'Op-Amp è spesso configurato come buffer (guadagno unitario).
        # Se si vuole un guadagno non unitario (tipico per filtri a campana o shelving con Op-Amp),
        # si usa un partitore resistivo in retroazione.
        # Per un Sallen-Key passa-basso, il guadagno è di solito unitario per un buon Q.
        # Se opamp_gain_resistor_ratio != 1.0, possiamo implementare un guadagno non unitario.
        # Per il momento, manteniamo la configurazione a buffer (V- connesso a Vout).
        # self.connect_nodes("OpAmp_Vminus", "OpAmp_Vout") # Questo configura il buffer

        # Se volessimo un guadagno A > 1:
        # A = 1 + R_f / R_g
        # dove R_f e R_g sono resistori di feedback sull'ingresso invertente.
        # Se self.opamp_gain_res_ratio == 1.0, significa A=1 (buffer).
        # Se self.opamp_gain_res_ratio > 1.0, usiamo un partitore.
        if self.opamp_gain_res_ratio > 1.0 + 1e-9: # Se il guadagno richiesto non è unitario
            # R_f e R_g sono calcolati per il guadagno A = 1 + R_f/R_g
            # Per semplicità, scegliamo R_g e calcoliamo R_f.
            R_g_val = 10000.0 # Resistenza al GND
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
        Formule per Sallen-Key passa-basso con guadagno unitario (A=1):
        Fc = 1 / (2 * pi * sqrt(R1*R2*C1*C2))
        Q = sqrt(R1*R2*C1*C2) / (C2 * (R1 + R2))
        Se R1=R2=R e C1=C2=C:
        Fc = 1 / (2 * pi * R * C)
        Q = 1 / (2) (per A=1, R1=R2, C1=C2) -> Q=0.5
        Se C1 = 2*C, C2 = C, R1=R2=R:
        Fc = 1 / (2 * pi * R * sqrt(2) * C)
        Q = 1 / sqrt(2) approx 0.707 (Butterworth)
        Per un guadagno A: Q = sqrt(R1*R2*C1*C2) / (C2 * (R1 + R2) + C1*(R1+R2)*(1-A)) (complesso)
        """
        # Assumiamo R1=R2=R e C1=C2=C per calcoli più semplici, per ora.
        # Per il modello attuale, usiamo le R e C individuali.

        # Frequenza di taglio
        if self.R1_val * self.R2_val * self.C1_val * self.C2_val <= 0:
            fc = 0.0 # Evita divisione per zero
        else:
            fc = 1.0 / (2.0 * np.pi * np.sqrt(self.R1_val * self.R2_val * self.C1_val * self.C2_val))

        # Fattore Q (per guadagno unitario)
        # Se il guadagno non è unitario, la formula del Q diventa più complessa e dipende anche da A.
        # Q = sqrt(R1*R2*C1*C2) / (C1*(R1+R2) + C2*R2*(1-A_v))
        # Dove A_v è il guadagno dell'Op-Amp (self.opamp_gain_res_ratio)
        A_v = self.opamp_gain_res_ratio
        
        denominator = (self.C1_val * (self.R1_val + self.R2_val)) + (self.C2_val * self.R1_val * (1 - A_v))
        
        if denominator == 0:
            Q = float('inf') # Risonanza infinita
        else:
            numerator = np.sqrt(self.R1_val * self.R2_val * self.C1_val * self.C2_val)
            Q = numerator / denominator

        return {"cutoff_frequency": fc, "Q": Q, "gain": A_v}

    def set_filter_parameters(self, cutoff_freq, Q_factor, gain=1.0):
        """
        Aggiorna i valori dei componenti R e C per raggiungere una data Fc e Q.
        Questo è il metodo "chiave" per rendere il filtro regolabile.
        Richiede di risolvere un sistema di equazioni.

        Per un Sallen-Key passa-basso:
        Assumiamo C1 = C2 = C_fixed
        Fc = 1 / (2 * pi * R * C) => R = 1 / (2 * pi * Fc * C)
        Q = 1 / (3 - A) => A = 3 - 1/Q (per un caso specifico R1=R2, C1=C2)

        Per semplicità, potremmo fissare R1, R2 e calcolare C1, C2.
        Oppure, fissare i condensatori e calcolare le resistenze.
        Dato che abbiamo 4 variabili (R1, R2, C1, C2) e 2 obiettivi (Fc, Q), abbiamo gradi di libertà.

        Una strategia comune:
        1. Scegli C1 e C2 (es. C1 = 1uF, C2 = 1uF)
        2. Calcola R1 e R2:
           R = 1 / (2 * pi * Fc * C)
           Per ottenere il Q desiderato, si può variare R1/R2 o C1/C2 o il guadagno A.
           Per un filtro di Sallen-Key, Q = 1 / (2 - A) (se R1C1 = R2C2)
           Questo non è il Sallen-Key più flessibile per variare il Q.

        Per una variazione "dolce" vs "campanatura", si regola il Q.
        Un Q più basso (es. 0.5) = più dolce (Bessel).
        Un Q = 0.707 = Butterworth (risposta piatta in banda passante).
        Un Q > 0.707 = Campanatura (Chebyshev).

        Questo è il punto più complesso. Per rendere l'Op-Amp "regolabile",
        potremmo variare le resistenze e i condensatori, ma in pratica si varierebbe
        solo 1-2 valori con potenziometri.

        Per il momento, lasciamo questo metodo come placeholder o per debug.
        L'utente imposterà i valori R1, R2, C1, C2 direttamente nella `__init__`.
        Se vogliamo "regolabilità", dobbiamo aggiungere "variabili" al circuito.

        Questo metodo è da implementare se si vuole una regolazione diretta da Fc, Q, Gain.
        """
        print(f"DEBUG: Tentativo di impostare Fc={cutoff_freq}Hz, Q={Q_factor}, Gain={gain} (non implementato completamente per ricalcolo R,C)")
        # La logica per ricalcolare R, C, A per un dato Fc, Q, Gain è complessa
        # e richiede un algoritmo di sintesi del filtro.
        # Per ora, è meglio impostare R,C all'inizio e osservare i loro parametri.
        pass
