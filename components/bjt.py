# components/bjt.py

import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM # Importa costanti fisiche
from utils.helpers import numerical_jacobian # Importa per le derivate numeriche

class BJT:
    def __init__(self, Is=1e-14, Bf=200, Br=1.0, Vaf=80.0, Var=20.0, Ise=1e-14, Isc=1e-14, Ne=1.5, Nc=2.0):
        """
        Modello semplificato di Transistor Bipolare a Giunzione (BJT) NPN.
        Basato su un sottoinsieme dei parametri del modello Gummel-Poon (Ebers-Moll semplificato).
        Args:
            Is (float): Corrente di saturazione inversa di giunzione (Ampere).
            Bf (float): Guadagno di corrente in forward (beta_F).
            Br (float): Guadagno di corrente in reverse (beta_R).
            Vaf (float): Tensione di Early in forward (Volt).
            Var (float): Tensione di Early in reverse (Volt).
            Ise (float): Corrente di saturazione per emettitore non ideale.
            Isc (float): Corrente di saturazione per collettore non ideale.
            Ne (float): Coefficiente di emissione per emettitore.
            Nc (float): Coefficiente di emissione per collettore.
        """
        self.Is = float(Is)
        self.Bf = float(Bf) # Forward Beta (hFE)
        self.Br = float(Br) # Reverse Beta
        self.Vaf = float(Vaf) # Forward Early Voltage
        self.Var = float(Var) # Reverse Early Voltage
        self.Ise = float(Ise) # Non-ideal base-emitter saturation current
        self.Ne = float(Ne) # Non-ideal base-emitter emission coefficient
        self.Isc = float(Isc) # Non-ideal base-collector saturation current
        self.Nc = float(Nc) # Non-ideal base-collector emission coefficient

        self.Vt = VT_ROOM # Tensione termica

    def calculate_collector_current(self, v_be, v_ce):
        """
        Calcola la corrente di Collettore (Ic) del BJT NPN.
        Modello Ebers-Moll semplificato con effetto Early.
        Args:
            v_be (float): Tensione Base-Emettitore (Vbe).
            v_ce (float): Tensione Collettore-Emettitore (Vce).
        Returns:
            float: Corrente di Collettore (Ic) in Ampere.
        """
        # Corrente di Emettitore ideale (forward)
        exp_vbe = np.exp(v_be / self.Vt)
        I_f = self.Is * (exp_vbe - 1.0) # Corrente ideale di diodo per Vbe

        # Corrente di Collettore ideale (reverse) - non sempre usata in modelli semplici
        # exp_vbc = np.exp((v_be - v_ce) / self.Vt)
        # I_r = self.Is * (exp_vbc - 1.0)

        # Effetto Early (modulazione di base)
        # 1 + Vce / Vaf (forward active)
        # 1 + Vbc / Var (reverse active)
        # Per forward active:
        early_factor_f = (1.0 + v_ce / self.Vaf) if self.Vaf > 0 else 1.0
        early_factor_r = (1.0 + (v_ce - v_be) / self.Var) if self.Var > 0 else 1.0 # Vbc = Vbe - Vce

        # Corrente di Collettore (attiva diretta)
        Ic = (self.Bf / (self.Bf + 1.0)) * I_f * early_factor_f

        # Aggiungi corrente di saturazione/cut-off
        if v_be < 0.5: # Sotto la tensione di accensione tipica
            Ic = 0.0 # Cut-off
        elif v_ce < 0.2 and Ic > 0: # Regione di saturazione (Vce basso)
            # In saturazione, Ic non segue il modello standard beta_f * Ib
            # Si avvicina a Vce_sat/R_sat. Qui si può fare una limitazione.
            # Questo è molto semplificato. Un modello SPICE avrebbe più equazioni.
            Ic = min(Ic, (v_ce / 10.0)) # Corrente limitata da una resistenza bassa (10 Ohm)
            pass

        return max(0.0, Ic) # La corrente non può essere negativa

    def calculate_base_current(self, v_be, v_ce):
        """
        Calcola la corrente di Base (Ib) del BJT NPN.
        """
        # Corrente di Emettitore ideale (forward)
        exp_vbe = np.exp(v_be / self.Vt)
        I_f = self.Is * (exp_vbe - 1.0)

        # Corrente di Base ideale (da corrente di Collettore)
        # Ib = Ic / Bf
        Ib = self.calculate_collector_current(v_be, v_ce) / self.Bf

        # Aggiungi corrente di ricombinazione (non ideale)
        # Ibe_recomb = Ise * (exp(Vbe / (Ne*Vt)) - 1)
        # Ib += Ibe_recomb / self.Bf # Semplificazione

        return max(0.0, Ib)

    # --- Metodi per la Jacobiana ---
    # Questi calcolano le derivate parziali necessarie per la matrice Jacobiana.
    # Per semplicità, qui usiamo la derivata numerica. Per precisione, implementare analiticamente.

    def calculate_transconductance(self, v_be, v_ce):
        """
        Calcola la transconduttanza (gm = d(Ic)/d(Vbe)) del BJT.
        """
        return numerical_jacobian(lambda v: self.calculate_collector_current(v, v_ce), v_be)

    def calculate_input_conductance(self, v_be, v_ce):
        """
        Calcola la conduttanza di ingresso (gpi = d(Ib)/d(Vbe)) del BJT.
        """
        return numerical_jacobian(lambda v: self.calculate_base_current(v, v_ce), v_be)

    def calculate_output_conductance(self, v_be, v_ce):
        """
        Calcola la conduttanza di uscita (gce = d(Ic)/d(Vce)) del BJT.
        """
        return numerical_jacobian(lambda v: self.calculate_collector_current(v_be, v), v_ce)

    def __str__(self):
        return f"BJT(Is={self.Is:.1e}, Bf={self.Bf:.0f}, Vaf={self.Vaf:.1f}V)"
