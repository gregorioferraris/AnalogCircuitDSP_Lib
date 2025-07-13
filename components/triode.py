# components/triode.py

import numpy as np
from utils.helpers import numerical_jacobian # Importa per le derivate numeriche

class Triode:
    def __init__(self, mu=100.0, Kp=600.0, X=1.5, Kg1=2.5):
        """
        Modello semplificato di Triodo a vuoto (es. based on Koren model/SPICE triode models).
        Args:
            mu (float): Fattore di amplificazione (gain).
            Kp (float): Parametro di perveanza (mA/V^(3/2) o mA/V^X).
                       A volte chiamato k.
            X (float): Esponente nel modello di corrente (solitamente 1.5 per le valvole).
            Kg1 (float): Parametro di controllo della griglia.
                         Spesso modellato come 1/mu + 1/gain_factor.
        """
        self.mu = float(mu) # Amplification factor
        self.Kp = float(Kp) # Perveance parameter (mA/V^X)
        self.X = float(X)   # Exponent (typically 1.5 for triodes)
        self.Kg1 = float(Kg1) # Grid control parameter (similar to Vgk / mu)

        # Converti Kp da mA/V^X a A/V^X se necessario (se la corrente è in Ampere)
        self.Kp_A = self.Kp / 1000.0 # Se Kp è in mA/V^X, converti ad A/V^X

    def calculate_anode_current(self, v_grid, v_anode):
        """
        Calcola la corrente di Anodo (Ia) del Triodo.
        Modello semplificato (es. Koren-like) per Ia.
        Args:
            v_grid (float): Tensione Griglia-Catodo (Vg).
            v_anode (float): Tensione Anodo-Catodo (Va).
        Returns:
            float: Corrente di Anodo (Ia) in Ampere.
        """
        # Tensione effettiva che controlla la corrente
        # Veff = (Va/mu + Vg)
        # Un altro modello: Veff = (Vg + Va/mu) / (1 + Va/(mu*Kg1))
        # O modello più semplice: Veff = (Vg + Va/mu)
        # La corrente è proporzionale a (Veff)^X

        # Un modello comune (es. per Simulink o Spice semplificato)
        # If Vg_eff > 0, then Ia = Kp * (Vg_eff)^X
        # Vg_eff = Vg + Va/mu (o con termini aggiuntivi per limitare la griglia)

        # Per il modello Koren:
        # V_grid_eff = (V_grid + V_anode / self.mu) / (1.0 + (V_anode / (self.mu * self.Kg1)))
        # Ia = self.Kp_A * (V_grid_eff**self.X)
        # Questo modello è complesso per la derivazione.

        # Usiamo una versione più semplice, ma che catturi la non linearità:
        # Ia = K * (V_g + V_a / mu)^X
        # Aggiungiamo una soglia di conduzione V_th per Vg+Va/mu
        V_eff = v_grid + v_anode / self.mu
        if V_eff <= 0: # Sotto la cut-off
            return 0.0
        
        # Gestione della corrente di griglia (solo se V_grid > 0, altrimenti è 0)
        # Questo non è la corrente di anodo, ma per modelli più complessi è importante.
        # Per ora ci concentriamo solo su Ia.
        
        # Aggiungi un piccolo termine per evitare problemi numerici vicino a 0
        ia_current = self.Kp_A * (V_eff**self.X) * (1 + 0.001 * v_anode) # Piccola dipendenza da Va anche in saturazione

        return max(0.0, ia_current) # Corrente di anodo non può essere negativa

    # --- Metodi per la Jacobiana ---
    # Utilizzano numerical_jacobian per semplicità, ma per performance e precisione
    # le derivate analitiche sono migliori.

    def calculate_transconductance(self, v_grid, v_anode):
        """
        Calcola la transconduttanza (gm = d(Ia)/d(Vg)) del Triodo.
        """
        return numerical_jacobian(lambda v: self.calculate_anode_current(v, v_anode), v_grid)

    def calculate_plate_conductance(self, v_grid, v_anode):
        """
        Calcola la conduttanza di placca (gp = d(Ia)/d(Va)) del Triodo.
        """
        return numerical_jacobian(lambda v: self.calculate_anode_current(v_grid, v), v_anode)

    def __str__(self):
        return f"Triode(mu={self.mu:.1f}, Kp={self.Kp:.1f} mA/V^{self.X:.1f})"
