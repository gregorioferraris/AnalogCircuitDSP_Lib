# components/mosfet.py

import numpy as np
from utils.helpers import numerical_jacobian # Importa per le derivate numeriche

class MOSFET:
    def __init__(self, name="MOSFET", Vt=1.0, Kn=1.0e-3, lambda_val=0.0):
        """
        Modello semplificato di MOSFET N-channel (può simulare enrichment o depletion-mode).
        Per simulare un JFET a impoverimento, impostare Vt con un valore negativo (e.g., -2.0 a -6.0V).
        Args:
            name (str): Nome identificativo del componente.
            Vt (float): Tensione di soglia (Threshold Voltage) (Volt).
                        Per N-channel Enhancement: Vt > 0.
                        Per N-channel Depletion (o JFET simile): Vt < 0.
            Kn (float): Parametro di transconduttanza (spesso Kp o Beta) (A/V^2).
                        Include (mu_n * Cox / 2) * (W/L).
            lambda_val (float): Parametro di modulazione di canale (per Id vs Vds in saturazione).
        """
        self.name = name
        self.Vt = float(Vt)
        self.Kn = float(Kn)
        self.lambda_val = float(lambda_val)

        # Mappatura dei nodi per l'integrazione nel circuito
        self.nodes = {} # Chiavi: 'drain', 'gate', 'source'

    def set_nodes(self, drain_node_id, gate_node_id, source_node_id):
        """Assegna gli ID dei nodi ai pin del MOSFET."""
        self.nodes['drain'] = drain_node_id
        self.nodes['gate'] = gate_node_id
        self.nodes['source'] = source_node_id

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
    # Per semplicità, qui usiamo la derivata numerica.

    def calculate_jacobian_elements(self, v_source, v_gate, v_drain):
        """
        Calcola gli elementi della matrice Jacobiana per il MOSFET.
        Restituisce una sottomatrice Jacobiana 3x3 per (Drain, Gate, Source)
        per la gestione nel solutore MNA.

        J_ij = d(I_i) / d(V_j)
        Dove I_i è la corrente che *esce* dal nodo i (standard MNA),
        e V_j è la tensione al nodo j.
        Le correnti I_Drain, I_Gate, I_Source.
        Per un MOSFET ideale, I_Gate = 0.
        I_Source = -(I_Drain + I_Gate) = -I_Drain.

        Quindi dobbiamo calcolare:
        d(Id)/dVd, d(Id)/dVg, d(Id)/dVs
        d(Ig)/dVd, d(Ig)/dVg, d(Ig)/dVs (= 0 per Ig=0)
        d(Is)/dVd, d(Is)/dVg, d(Is)/dVs
        """
        # Tensioni locali
        v_gs = v_gate - v_source
        v_ds = v_drain - v_source

        # Matrice Jacobiana 3x3 per i nodi (Drain, Gate, Source)
        jacobian = np.zeros((3, 3))

        # Funzione helper per calcolare le derivate numeriche
        # f(x_val, y_val) = self.calculate_drain_current(v_gs=x_val, v_ds=y_val)
        # Per la Jacobiana, le derivate sono rispetto alle tensioni dei nodi.
        
        # d(Id)/d(Vd) = d(Id)/d(Vds) * d(Vds)/d(Vd) = gds * 1 = gds
        gds = numerical_jacobian(lambda v: self.calculate_drain_current(v_gs, v), v_ds)
        
        # d(Id)/d(Vg) = d(Id)/d(Vgs) * d(Vgs)/d(Vg) = gm * 1 = gm
        gm = numerical_jacobian(lambda v: self.calculate_drain_current(v, v_ds), v_gs)

        # d(Id)/d(Vs) = d(Id)/d(Vgs) * d(Vgs)/d(Vs) + d(Id)/d(Vds) * d(Vds)/d(Vs)
        #             = gm * (-1) + gds * (-1) = -(gm + gds)
        
        # Riempimento della matrice Jacobiana
        # Righe: Correnti dai nodi (I_Drain, I_Gate, I_Source)
        # Colonne: Variazione delle tensioni ai nodi (V_Drain, V_Gate, V_Source)

        # Riga 0: Corrente di Drain (Id)
        jacobian[0, 0] = gds           # d(Id)/d(Vd)
        jacobian[0, 1] = gm            # d(Id)/d(Vg)
        jacobian[0, 2] = -(gm + gds)   # d(Id)/d(Vs)

        # Riga 1: Corrente di Gate (Ig) - Assumiamo Ig = 0 per MOSFET ideale
        jacobian[1, 0] = 0.0
        jacobian[1, 1] = 0.0
        jacobian[1, 2] = 0.0

        # Riga 2: Corrente di Source (Is = -Id)
        jacobian[2, 0] = -gds          # d(-Id)/d(Vd)
        jacobian[2, 1] = -gm           # d(-Id)/d(Vg)
        jacobian[2, 2] = (gm + gds)    # d(-Id)/d(Vs)

        return jacobian

    def __str__(self):
        return f"MOSFET({self.name}, Vt={self.Vt:.1f}V, Kn={self.Kn:.1e}, lambda={self.lambda_val:.2f})"
