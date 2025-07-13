# circuits/MultipleFeedbackBandPassCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.op_amp import OpAmp

import numpy as np

class MultipleFeedbackBandPassCircuit(Circuit):
    """
    Circuito di un filtro passa-banda attivo di tipo Multiple Feedback (MFB).
    Questo filtro è spesso usato come filtro a campana (Peak/Bell) in equalizzatori parametrici,
    permettendo di controllare la frequenza centrale, il Q e il guadagno.
    È un filtro del secondo ordine.
    """
    def __init__(self, name="MFB_BandPass_Bell", R1=10000.0, R2=10000.0, R3=10000.0,
                 C1=1.0e-8, C2=1.0e-8, opamp_params=None, sample_rate=48000):
        super().__init__(name)
        self.R1_val = float(R1)
        self.R2_val = float(R2)
        self.R3_val = float(R3) # Resistenza di guadagno (per il picco)
        self.C1_val = float(C1)
        self.C2_val = float(C2)
        self.opamp_params = opamp_params if opamp_params is not None else {}
        self.sample_rate = sample_rate

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il filtro MFB."""
        self.add_node("Input")
        self.add_node("Output")
        self.add_node("OpAmp_Vout")
        self.add_node("OpAmp_Vminus") # Ingresso invertente (sommatore)
        self.add_node("OpAmp_Vplus")  # Ingresso non invertente (a GND)

        # Nodi interni
        self.add_node("Node_A") # Nodo tra C1, R2, R3
        self.add_node("Node_B") # Nodo tra R1, C1

    def _add_components(self):
        """Aggiunge i componenti (R, C, Op-Amp)."""
        # Resistor R1 (ingresso)
        self.add_component(Resistor(self.R1_val), "Input", "OpAmp_Vminus")

        # Capacitor C1 (feedback)
        self.add_component(Capacitor(self.C1_val, sample_rate=self.sample_rate), "OpAmp_Vminus", "Node_A")

        # Resistor R2 (feedback)
        self.add_component(Resistor(self.R2_val), "Node_A", "OpAmp_Vout")

        # Capacitor C2 (feedback)
        self.add_component(Capacitor(self.C2_val, sample_rate=self.sample_rate), "Node_A", "GND")

        # Resistor R3 (per il guadagno, dal GND all'ingresso invertente)
        self.add_component(Resistor(self.R3_val), "OpAmp_Vminus", "GND")


        # Op-Amp
        # L'ingresso non invertente è a GND virtuale per questa topologia.
        self.op_amp = OpAmp(**self.opamp_params)
        self.add_component(self.op_amp, "OpAmp_Vplus", "OpAmp_Vminus", "OpAmp_Vout")

    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        # Ingresso non invertente dell'Op-Amp a GND (virtual ground per questo circuito)
        self.connect_nodes("OpAmp_Vplus", "GND")
        # L'output del filtro è l'uscita dell'Op-Amp
        self.connect_nodes("Output", "OpAmp_Vout")

    def calculate_parameters(self):
        """
        Calcola la frequenza centrale (Fc), il fattore di merito (Q) e il guadagno del picco.
        Formule per filtro MFB Band-Pass:
        fc = 1 / (2 * pi * sqrt(R2 * R3 * C1 * C2))
        Q = 1 / (2 * pi * fc * C1 * (R1 * R2 * R3 / (R1*R2 + R2*R3 + R1*R3)))
        Guadagno (a Fc): -R2 / R1
        Assumendo C1=C2=C per semplificare la sintesi per il controllo del Q.
        Se C1=C2=C: fc = 1 / (2 * pi * C * sqrt(R2 * R3))
                   Q = sqrt(R2 * R3) / (R1 * C * (1 + R2/R1 + R2/R3)) # semplificando
        Se R1, R2, R3 sono uguali e C1, C2 uguali, il Q non è facilmente regolabile
        per un dato Fc e Gain. Per la flessibilità, R1, R2, R3, C1, C2 devono essere scelti attentamente.
        """
        # Per semplicità, usiamo le formule per C1=C2=C
        C_eff = np.sqrt(self.C1_val * self.C2_val) # Media geometrica se diverse
        if self.R2_val * self.R3_val <= 0 or C_eff <= 0:
            fc = 0.0
        else:
            fc = 1.0 / (2.0 * np.pi * C_eff * np.sqrt(self.R2_val * self.R3_val))

        gain_at_fc_linear = -self.R2_val / self.R1_val # Guadagno è negativo per questa topologia
        gain_at_fc_db = 20 * np.log10(abs(gain_at_fc_linear))

        # Calcolo del Q: più complesso e dipende da tutti i valori.
        # Una formula semplificata per Q in MFB se R3 (ingresso non invertente) è a GND
        # Q = 1 / (2 * pi * fc * C1 * R1 * (1 + R2/R1)) # approssimazione
        
        # Formula più generale per MFB BandPass:
        # A_v = -R2 / R1
        # Q = 1 / (sqrt(C1/C2) * ( (R1+R2)/(R1*R2) * sqrt(R2*R3) ) + 1/R3 * sqrt(R2*R3*C1/C2) )
        # Questo è un po' più complesso da generalizzare.
        # In pratica, si fissa C1=C2, si sceglie R1 per il guadagno, e R2/R3 per Fc e Q.
        # R2 = Q / (2 * pi * fc * C)
        # R3 = Q / (2 * pi * fc * C * (Q/Gain - 1))
        # R1 = R2 / |Gain|
        
        # Per ora, restituiamo un placeholder per Q o una formula semplificata
        # Per una sintesi diretta da Fc, Q, Gain, l'utente dovrebbe fornire C1 e C2 fissi.
        # Esempio per Q, assumendo C1=C2=C e R1, R2, R3 come calcolati per Fc, Q, Gain:
        # Questo richiede una funzione di sintesi inversa.
        # Per questa implementazione diretta, daremo un Q approssimato o nullo.
        
        # Per un Q accurato, le formule sono complesse e di solito si assumono C1=C2=C
        # e si calcolano le R per ottenere Fc, Q, Gain.
        # Q_val = (1.0 / (2.0 * np.pi * fc * self.C1_val)) * np.sqrt(1.0/(self.R1_val) + 1.0/self.R2_val + 1.0/self.R3_val)
        # Il Q per MFB è: sqrt(R_eff * C_eff) / (R1*C1 + R2*C2)
        # Con R1, R2, R3, C1, C2 generali, il calcolo del Q è difficile senza assumere relazioni.
        # Userò una stima approssimativa o lascerò come da calcolare esternamente.
        
        # Visto che R2 e R3 sono nel loop di feedback e C1, C2 sono cruciali.
        # Per avere un controllo "dolce" vs "campanatura", si agisce tipicamente sul Q.
        # Un modo per sintetizzare R, C per MFB:
        # Se C1 = C2 = C:
        # K = 2 * pi * fc * C
        # R1 = abs(Gain) / K
        # R2 = 1 / (K * (1/Q - Gain))  (se Gain negativo)
        # R3 = 1 / (K * (1/Q))
        
        # Poiché l'utente imposterà R1, R2, R3, C1, C2, calcoleremo il Q risultante
        # in base a queste scelte.
        # Per MFB, se R1=R, R2=R, R3=R e C1=C2=C:
        # Fc = 1 / (2*pi*R*C)
        # Q = 1 / 3
        # Per un Q maggiore (campanatura), si gioca con R2 e R3 in relazione a R1.
        
        # Riportiamo semplicemente i valori nominali, il Q sarà intrinseco.
        return {"center_frequency": fc, "Q": "Determined by R, C values", "gain_db": gain_at_fc_db}

    def set_parameters(self, center_freq, Q_factor, gain_db=0.0):
        """
        Metodo per sintetizzare i valori di R e C per raggiungere Fc, Q, e Gain.
        Richiede un approccio di sintesi e non è un semplice aggiornamento dei componenti.
        Questo è un placeholder. In un sistema reale, potresti avere potenziometri
        che variano R o C, o una logica che calcola i valori dei componenti.

        Per la sintesi, tipicamente si fissa C1=C2=C e si calcolano R1, R2, R3.
        """
        print(f"DEBUG: Sintesi di MFB BandPass per Fc={center_freq}Hz, Q={Q_factor}, Gain={gain_db}dB (non implementato)")
        # Esempio di sintesi per MFB (semplificato, assumendo C1=C2=C_ref)
        # C_ref = 1.0e-8 # Fissa un valore di capacità di riferimento
        # K = 2 * np.pi * center_freq * C_ref
        # R1 = abs(10**(gain_db/20.0)) / K
        # R2 = 1.0 / (K * (1.0/Q_factor - 10**(gain_db/20.0)))
        # R3 = 1.0 / (K * (1.0/Q_factor))
        # Aggiornare self.R1_val, self.R2_val, self.R3_val, self.C1_val, self.C2_val.
        pass
