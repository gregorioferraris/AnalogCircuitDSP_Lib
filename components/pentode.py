# components/pentode.py

import numpy as np
from utils.helpers import numerical_jacobian # Importa per le derivate numeriche

class Pentode:
    def __init__(self, name="Pentode", mu=10.0, Kp=300.0, X=1.5, Kg1=5.0, Kg2=10.0):
        """
        Modello semplificato di Pentodo a vuoto (es. basato su modelli empirici).
        Args:
            name (str): Nome del componente.
            mu (float): Fattore di amplificazione (mu della triode-section implicita).
            Kp (float): Parametro di perveanza.
            X (float): Esponente nel modello di corrente (solitamente 1.5).
            Kg1 (float): Parametro di controllo della griglia di controllo (g1).
            Kg2 (float): Parametro di controllo della griglia schermo (g2).
        """
        self.name = name # Aggiunto il nome
        self.mu = float(mu)
        self.Kp = float(Kp) # Perveance parameter (mA/V^X)
        self.X = float(X)
        self.Kg1 = float(Kg1) # Grid 1 control factor
        self.Kg2 = float(Kg2) # Grid 2 control factor

        self.Kp_A = self.Kp / 1000.0 # Converti Kp da mA/V^X a A/V^X

        # Mappatura dei nodi (Anodo/Placca, Griglia1, Griglia2, Catodo)
        self.nodes = {} # Chiavi: 'anode', 'grid1', 'grid2', 'cathode'

    def set_nodes(self, anode_node_id, grid1_node_id, grid2_node_id, cathode_node_id):
        """Assegna gli ID dei nodi ai pin del Pentodo."""
        self.nodes['anode'] = anode_node_id
        self.nodes['grid1'] = grid1_node_id
        self.nodes['grid2'] = grid2_node_id
        self.nodes['cathode'] = cathode_node_id

    def calculate_anode_current(self, v_grid1_cathode, v_grid2_cathode, v_anode_cathode):
        """
        Calcola la corrente di Anodo (Ia) del Pentodo.
        Modello semplificato che considera l'influenza di g1, g2 e anodo.
        Args:
            v_grid1_cathode (float): Tensione Griglia-1-Catodo (Vg1k).
            v_grid2_cathode (float): Tensione Griglia-2-Catodo (Vg2k).
            v_anode_cathode (float): Tensione Anodo-Catodo (Vak).
        Returns:
            float: Corrente di Anodo (Ia) in Ampere.
        """
        # Modello basato su Vg1 e Vg2. Questo modello assume che il catodo sia il riferimento (GND).
        # Le tensioni sono relative al catodo.

        # Tensione efficace di griglia (combinazione di Vg1 e Vg2)
        # Questo è il punto critico del modello.
        # Ho adattato la tua V_eff dalla tua implementazione originale:
        V_eff_g = v_grid1_cathode + v_grid2_cathode / self.Kg2 # Tua V_eff

        # Se V_eff_g è negativa o zero, la valvola è in cut-off.
        # (Questo è un cut-off semplificato, i pentodi hanno curve più complesse)
        if V_eff_g <= 0:
            return 0.0

        # Termine di modulazione di canale (plate resistance)
        # Per pentodi in regione di saturazione, la corrente è relativamente indipendente da Va.
        # Tuttavia, c'è sempre una leggera pendenza (resistenza di placca finita).
        # Questo è il lambda dalla tua implementazione.
        lambda_anode = 0.01 # Un valore molto piccolo, come nel tuo codice.

        ia_current = self.Kp_A * (V_eff_g**self.X) * (1.0 + lambda_anode * v_anode_cathode)
        
        return max(0.0, ia_current) # Assicurati che la corrente non sia mai negativa

    # --- Metodi per la Jacobiana ---
    # Utilizzano numerical_jacobian per semplicità e generalità.

    def calculate_jacobian_elements(self, v_anode, v_grid1, v_grid2, v_cathode):
        """
        Calcola gli elementi della matrice Jacobiana per il Pentodo.
        Restituisce una sottomatrice Jacobiana 4x4 per (Anodo, Griglia1, Griglia2, Catodo).

        J_ij = d(I_i) / d(V_j)
        Dove I_i è la corrente che *esce* dal nodo i (standard MNA),
        e V_j è la tensione al nodo j.
        Correnti: I_Anodo, I_Griglia1, I_Griglia2, I_Catodo.
        I_Griglia1 e I_Griglia2 sono solitamente molto piccole o zero per i pentodi,
        a meno di modellare la corrente di griglia. Per ora, le assumiamo zero.
        I_Catodo = -(I_Anodo + I_Griglia1 + I_Griglia2) = -I_Anodo.
        """
        # Tensioni locali rispetto al catodo
        v_grid1_cathode = v_grid1 - v_cathode
        v_grid2_cathode = v_grid2 - v_cathode
        v_anode_cathode = v_anode - v_cathode

        # Matrice Jacobiana 4x4 per i nodi (Anodo, Griglia1, Griglia2, Catodo)
        jacobian = np.zeros((4, 4))

        # Calcolo delle derivate parziali necessarie usando numerical_jacobian
        # gm1 = d(Ia)/d(Vg1k)
        gm1 = numerical_jacobian(lambda v: self.calculate_anode_current(v, v_grid2_cathode, v_anode_cathode), v_grid1_cathode)
        
        # gm2 = d(Ia)/d(Vg2k)
        gm2 = numerical_jacobian(lambda v: self.calculate_anode_current(v_grid1_cathode, v, v_anode_cathode), v_grid2_cathode)
        
        # gp = d(Ia)/d(Vak)
        gp = numerical_jacobian(lambda v: self.calculate_anode_current(v_grid1_cathode, v_grid2_cathode, v), v_anode_cathode)

        # Riempimento della matrice Jacobiana
        # Righe: Correnti dai nodi (I_Anode, I_Grid1, I_Grid2, I_Cathode)
        # Colonne: Variazione delle tensioni ai nodi (V_Anode, V_Grid1, V_Grid2, V_Cathode)

        # Riga 0: Corrente di Anodo (Ia)
        jacobian[0, 0] = gp                     # d(Ia)/d(Va)
        jacobian[0, 1] = gm1                    # d(Ia)/d(Vg1)
        jacobian[0, 2] = gm2                    # d(Ia)/d(Vg2)
        jacobian[0, 3] = -(gm1 + gm2 + gp)      # d(Ia)/d(Vk) (dal KCL Ia + Ig1 + Ig2 + Ik = 0)

        # Righe 1 e 2: Correnti di Griglia (Ig1, Ig2) - Assumiamo Ig1 = 0, Ig2 = 0 per ora
        # Se volessimo modellare Ig1 e Ig2 (per Vg1k > 0 o Vg2k > 0), avremmo bisogno di
        # funzioni calculate_grid1_current() e calculate_grid2_current() e le loro Jacobian.
        jacobian[1, :] = 0.0 # Riga per I_Grid1
        jacobian[2, :] = 0.0 # Riga per I_Grid2

        # Riga 3: Corrente di Catodo (Ik = -(Ia + Ig1 + Ig2) = -Ia per Ig=0)
        jacobian[3, 0] = -gp                    # d(-Ia)/d(Va)
        jacobian[3, 1] = -gm1                   # d(-Ia)/d(Vg1)
        jacobian[3, 2] = -gm2                   # d(-Ia)/d(Vg2)
        jacobian[3, 3] = (gm1 + gm2 + gp)       # d(-Ia)/d(Vk)

        return jacobian

    def __str__(self):
        return f"Pentode({self.name}, mu={self.mu:.1f}, Kp={self.Kp:.1f} mA/V^{self.X:.1f})"
