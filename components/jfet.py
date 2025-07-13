# components/jfet.py

import numpy as np
from utils.helpers import numerical_jacobian # Importa per le derivate numeriche

class JFET:
    def __init__(self, Idss=0.01, Vp=-2.0, lambda_val=0.05):
        """
        Modello semplificato di JFET (Junction Field-Effect Transistor) N-channel.
        Args:
            Idss (float): Corrente di Drain-Source a Vgs=0 e Vds > |Vp| (Ampere).
            Vp (float): Tensione di Pinch-Off (Volt). Per N-channel è negativo.
            lambda_val (float): Parametro di modulazione di canale (per Id vs Vds in saturazione).
        """
        if Idss <= 0:
            raise ValueError("Idss deve essere positivo.")
        if Vp >= 0: # Per JFET N-channel, Vp deve essere negativo
            raise ValueError("Vp (tensione di pinch-off) deve essere negativo per JFET N-channel.")
        if lambda_val < 0:
            raise ValueError("Lambda (modulazione di canale) deve essere non negativo.")

        self.Idss = float(Idss)
        self.Vp = float(Vp) # Tensione di pinch-off
        self.lambda_val = float(lambda_val) # Parametro lambda

    def calculate_drain_current(self, v_gs, v_ds):
        """
        Calcola la corrente di Drain (Id) del JFET.
        Utilizza un modello semplificato che considera regione di triodo/saturazione.
        Args:
            v_gs (float): Tensione Gate-Source (Vgs).
            v_ds (float): Tensione Drain-Source (Vds).
        Returns:
            float: Corrente di Drain (Id) in Ampere.
        """
        if v_gs >= 0:
            # JFETs in genere non operano con Vgs positive (o sono limitate a pochi mV).
            # Per semplicità, consideriamo corrente nulla o molto bassa.
            # In realtà, ci sarebbe conduzione della giunzione gate-source.
            return 0.0 # Or some large current due to forward bias of gate junction

        if v_gs <= self.Vp:
            # Regione di Cut-off (Vgs < Vp)
            return 0.0
        elif v_ds >= (v_gs - self.Vp):
            # Regione di Saturazione (Vds >= Vgs - Vp)
            # Id = Idss * (1 - Vgs/Vp)^2 * (1 + lambda * Vds)
            term_vgs = (1.0 - v_gs / self.Vp)
            # Clamping per evitare NaN o numeri complessi se term_vgs diventa negativo a causa di approssimazioni
            if term_vgs < 0: term_vgs = 0
            id_saturation = self.Idss * (term_vgs ** 2) * (1.0 + self.lambda_val * v_ds)
            return max(0.0, id_saturation) # La corrente non può essere negativa
        else:
            # Regione di Triodo (Ohmica) (Vds < Vgs - Vp)
            # Modello semplificato: Id = Idss * [2(1 - Vgs/Vp)Vds - Vds^2]/Vp^2
            # Questo modello può essere instabile. Una linearizzazione è preferibile.
            # Una forma alternativa comune è Id = Beta * [2(Vgs-Vt)Vds - Vds^2]
            # Usiamo un'approssimazione lineare per la regione ohmica se non vogliamo il modello completo
            # Per semplicità, approssimiamo con una resistenza equivalente se Vds è piccolo.
            # Questo è un modello semplificato. Un modello più robusto userebbe la formula quadrata.
            # Per il momento, facciamo una transizione più o meno liscia.
            V_prime = v_gs - self.Vp
            id_triode = self.Idss * (2 * (1 - v_gs / self.Vp) * v_ds - v_ds**2 / self.Vp) / self.Vp
            return max(0.0, id_triode) # La corrente non può essere negativa

    def calculate_transconductance(self, v_gs, v_ds):
        """
        Calcola la transconduttanza (gm = d(Id)/d(Vgs)) del JFET.
        """
        # Se possibile, implementare la derivata analitica.
        # Per un modello semplificato, usiamo la derivata numerica.
        return numerical_jacobian(lambda v: self.calculate_drain_current(v, v_ds), v_gs)

    def calculate_output_conductance(self, v_gs, v_ds):
        """
        Calcola la conduttanza di output (gds = d(Id)/d(Vds)) del JFET.
        """
        # Se possibile, implementare la derivata analitica.
        # Per un modello semplificato, usiamo la derivata numerica.
        return numerical_jacobian(lambda v: self.calculate_drain_current(v_gs, v), v_ds)

    def __str__(self):
        return f"JFET(Idss={self.Idss:.1e}, Vp={self.Vp:.1f}V, lambda={self.lambda_val:.2f})"
