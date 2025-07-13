# components/pentode.py

import numpy as np
from utils.helpers import numerical_jacobian # Importa per le derivate numeriche

class Pentode:
    def __init__(self, mu=10.0, Kp=300.0, X=1.5, Kg1=5.0, Kg2=10.0):
        """
        Modello semplificato di Pentodo a vuoto (es. basato su modelli empirici).
        Args:
            mu (float): Fattore di amplificazione (mu della triode-section implicita).
            Kp (float): Parametro di perveanza.
            X (float): Esponente nel modello di corrente (solitamente 1.5).
            Kg1 (float): Parametro di controllo della griglia di controllo (g1).
            Kg2 (float): Parametro di controllo della griglia schermo (g2).
        """
        self.mu = float(mu)
        self.Kp = float(Kp) # Perveance parameter (mA/V^X)
        self.X = float(X)
        self.Kg1 = float(Kg1) # Grid 1 control factor
        self.Kg2 = float(Kg2) # Grid 2 control factor

        self.Kp_A = self.Kp / 1000.0 # Converti Kp da mA/V^X a A/V^X

    def calculate_anode_current(self, v_grid1, v_grid2, v_anode):
        """
        Calcola la corrente di Anodo (Ia) del Pentodo.
        Modello semplificato che considera l'influenza di g1, g2 e anodo.
        Args:
            v_grid1 (float): Tensione Griglia-1-Catodo (Vg1).
            v_grid2 (float): Tensione Griglia-2-Catodo (Vg2).
            v_anode (float): Tensione Anodo-Catodo (Va).
        Returns:
            float: Corrente di Anodo (Ia) in Ampere.
        """
        # Modello di Koren/D.T. Smith (semplificato per Ia)
        # Questo è un modello comune per i pentodi, ma richiede un'implementazione attenta.
        # Ia = Kp * (V_eff)^X * (1 + lambda * Va)
        # Dove V_eff = (Vg1 + Vg2/K + Va/mu)
        # O V_eff = (Vg1/Kg1 + Vg2/Kg2) se si usano coefficienti per le griglie
        
        # Una forma comune che include schermo e controllo:
        # V_eff = (v_grid1 + v_grid2 / self.Kg2) # Semplificazione per V_eff
        # Ia = self.Kp_A * (V_eff)**self.X # Questo ignora l'anodo.
        
        # Modello che include l'effetto di Va (plate resistance) e screengrid (g2)
        # V_grid_equivalent = v_grid1 + v_grid2 / self.mu_screen # mu_screen è un rapporto Vg2/Vg1 per Ia
        # Usiamo un modello più tipico:
        
        # Effetto di cut-off
        if v_grid1 <= 0: # Vg1 sotto cutoff, corrente zero
             return 0.0

        # Tensione efficace di griglia (combinazione di Vg1 e Vg2)
        # Questo è un punto chiave e può variare a seconda del modello.
        # Un approccio è trattare Vg2 come una sorgente di tensione ausiliaria costante.
        # E_g = (Vg1 + Vg2 / mu_G2) dove mu_G2 è fattore di amplificazione da G2 a Catodo
        # Per semplicità usiamo un modello in cui Vg2 influenza Kp
        
        # Variazione di Kp con Vg2 (molto semplificato)
        # Kp_adjusted = self.Kp_A * (1 + v_grid2 / self.Kg2) # Aumenta Kp con Vg2

        # Approccio basato su curve caratteristiche:
        # Corrente in saturazione di pentodo è relativamente indipendente da Va,
        # ma dipende da Vg1 e Vg2.
        # Ia = f(Vg1, Vg2)
        # Poi effetto di modulazione di canale: (1 + lambda * Va)
        
        # Proviamo con un modello basato su Vg1 e Vg2 (approccio di D.T. Smith)
        # Assume Va alta
        
        # Se (Vg1 + Vg2/Kg2) è la V_eff che guida Ia
        V_eff = v_grid1 + v_grid2 / self.Kg2

        if V_eff <= 0:
            return 0.0

        # Termine di modulazione di canale (plate resistance)
        # Lambda per i pentodi è molto più piccolo che per i triodi
        # Per pentodi in regione di saturazione, la corrente è quasi indipendente da Va.
        # Aggiungiamo un piccolo coefficiente lambda per l'anodo.
        lambda_anode = 0.01 # Un valore molto piccolo

        ia_current = self.Kp_A * (V_eff**self.X) * (1.0 + lambda_anode * v_anode)
        
        # Gestione della regione di triodo (se Va è molto bassa)
        # Quando Va scende sotto un certo livello (es. Vg2), il pentodo può entrare in triodo-region
        # Questa è una semplificazione del modello, per un pentodo "ideale" è sempre in saturazione
        # rispetto all'anodo se Va è sufficientemente alta.
        # Inizialmente non modelliamo la regione di triodo del pentodo.

        return max(0.0, ia_current)

    # --- Metodi per la Jacobiana ---
    # Utilizzano numerical_jacobian per semplicità.

    def calculate_transconductance_g1(self, v_grid1, v_grid2, v_anode):
        """
        Calcola la transconduttanza rispetto a Grid 1 (gm1 = d(Ia)/d(Vg1)).
        """
        return numerical_jacobian(lambda v: self.calculate_anode_current(v, v_grid2, v_anode), v_grid1)

    def calculate_transconductance_g2(self, v_grid1, v_grid2, v_anode):
        """
        Calcola la transconduttanza rispetto a Grid 2 (gm2 = d(Ia)/d(Vg2)).
        """
        return numerical_jacobian(lambda v: self.calculate_anode_current(v_grid1, v, v_anode), v_grid2)

    def calculate_plate_conductance(self, v_grid1, v_grid2, v_anode):
        """
        Calcola la conduttanza di placca (gp = d(Ia)/d(Va)) del Pentodo.
        """
        return numerical_jacobian(lambda v: self.calculate_anode_current(v_grid1, v_grid2, v), v_anode)

    def __str__(self):
        return f"Pentode(mu={self.mu:.1f}, Kp={self.Kp:.1f} mA/V^{self.X:.1f})"
