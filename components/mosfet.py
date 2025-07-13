# components/mosfet.py

import numpy as np
from utils.helpers import numerical_jacobian # Importa per le derivate numeriche

class MOSFET:
    def __init__(self, Vt=1.0, Kn=1.0e-3, lambda_val=0.0):
        """
        Modello semplificato di MOSFET N-channel ad arricchimento (Enhancement-mode).
        Args:
            Vt (float): Tensione di soglia (Threshold Voltage) (Volt).
            Kn (float): Transconduttanza di processo (spesso Kp o Beta) (A/V^2).
                       Include (mu_n * Cox / 2) * (W/L).
            lambda_val (float): Parametro di modulazione di canale (per Id vs Vds in saturazione).
        """
        self.Vt = float(Vt) # Threshold voltage (Vth)
        self.Kn = float(Kn) # Transconductance parameter (beta/2 in some models)
        self.lambda_val = float(lambda_val) # Channel-length modulation parameter

    def calculate_drain_current(self, v_gs, v_ds):
        """
        Calcola la corrente di Drain (Id) del MOSFET N-channel.
        Modello a grandi segnali che considera regione di Cut-off, Triodo (Lineare) e Saturazione.
        Args:
            v_gs (float): Tensione Gate-Source (Vgs).
            v_ds (float): Tensione Drain-Source (Vds).
        Returns:
            float: Corrente di Drain (Id) in Ampere.
        """
        # Effettiva tensione Gate-Source relativa alla soglia
        V_ov = v_gs - self.Vt # Overdrive voltage

        if V_ov <= 0:
            # Regione di Cut-off (Vgs <= Vt)
            return 0.0
        elif v_ds >= V_ov:
            # Regione di Saturazione (Vds >= Vgs - Vt)
            # Id = Kn * (Vgs - Vt)^2 * (1 + lambda * Vds)
            id_saturation = self.Kn * (V_ov ** 2) * (1.0 + self.lambda_val * v_ds)
            return max(0.0, id_saturation) # La corrente non può essere negativa
        else:
            # Regione di Triodo (Lineare/Ohmica) (Vds < Vgs - Vt)
            # Id = Kn * [2 * (Vgs - Vt) * Vds - Vds^2] * (1 + lambda * Vds)
            id_triode = self.Kn * (2 * V_ov * v_ds - v_ds**2) * (1.0 + self.lambda_val * v_ds)
            return max(0.0, id_triode) # La corrente non può essere negativa

    # --- Metodi per la Jacobiana ---
    # Questi calcolano le derivate parziali necessarie per la matrice Jacobiana.
    # Per semplicità, qui usiamo la derivata numerica. Per precisione, implementare analiticamente.

    def calculate_transconductance(self, v_gs, v_ds):
        """
        Calcola la transconduttanza (gm = d(Id)/d(Vgs)) del MOSFET.
        """
        # Derivata analitica in saturazione: gm = 2 * Kn * (Vgs - Vt) * (1 + lambda * Vds)
        # Per la regione di triodo e cut-off, la derivata analitica è diversa.
        # Usiamo la derivata numerica per coprire tutte le regioni.
        return numerical_jacobian(lambda v: self.calculate_drain_current(v, v_ds), v_gs)

    def calculate_output_conductance(self, v_gs, v_ds):
        """
        Calcola la conduttanza di output (gds = d(Id)/d(Vds)) del MOSFET.
        """
        # Derivata analitica in saturazione: gds = lambda * Kn * (Vgs - Vt)^2
        # Per la regione di triodo, la derivata analitica è diversa.
        # Usiamo la derivata numerica per coprire tutte le regioni.
        return numerical_jacobian(lambda v: self.calculate_drain_current(v_gs, v), v_ds)

    def __str__(self):
        return f"MOSFET(Vt={self.Vt:.1f}V, Kn={self.Kn:.1e}, lambda={self.lambda_val:.2f})"
