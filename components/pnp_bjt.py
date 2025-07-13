# components/pnp_bjt.py

import numpy as np
from utils.helpers import numerical_jacobian # Assicurati che questo percorso sia corretto e la funzione esista

class PNP_BJT:
    """
    Modello semplificato di BJT PNP.
    Basato sul modello Ebers-Moll semplificato per transistor PNP.
    """
    def __init__(self, name: str, nodes: dict,
                 Is: float = 1.0e-14,    # A (Corrente di saturazione inversa del diodo BE)
                 Vt: float = 0.02585,    # V (Tensione termica a 300K)
                 Beta_F: float = 100.0,  # hFE (Guadagno di corrente in polarizzazione diretta)
                 Va: float = 70.0):      # V (Tensione di Early)

        if Is <= 0 or Vt <= 0 or Beta_F <= 0 or Va <= 0:
            raise ValueError("I parametri Is, Vt, Beta_F e Va del BJT PNP devono essere positivi.")

        self.name = name
        self.nodes = nodes
        self.Is = float(Is)
        self.Vt = float(Vt)
        self.Beta_F = float(Beta_F)
        self.Va = float(Va)

        required_nodes = ['collector', 'base', 'emitter']
        if not all(node in self.nodes for node in required_nodes):
            raise ValueError(f"PNP_BJT '{name}': Nodi richiesti (collector, base, emitter) non specificati correttamente.")

        # Inizializza lo stato per il calcolo della Jacobiana
        self.last_Vbe = 0.0
        self.last_Vce = 0.0

    def calculate_collector_current(self, v_be: float, v_ce: float) -> float:
        """
        Calcola la corrente di Collettore (Ic) del BJT PNP.
        La corrente Ic è negativa se entra nel collettore e positiva se esce.
        Per un PNP, la corrente di collettore *esce* dal collettore (Ic < 0).

        Args:
            v_be (float): Tensione Base-Emettitore (V_B - V_E). Deve essere negativa per la conduzione.
            v_ce (float): Tensione Collettore-Emettitore (V_C - V_E). Deve essere negativa in regione attiva/saturazione.
        Returns:
            float: Corrente di Collettore (Ic) in Ampere (valore negativo quando il PNP conduce).
        """
        self.last_Vbe = v_be
        self.last_Vce = v_ce

        # Condizioni per il PNP
        # Cut-off: Se V_BE >= -0.5V (approssimazione) o più rigorosamente V_BE >= 0V (se non c'è tensione di accensione)
        # Il diodo BE deve essere polarizzato inversamente (Vbe positiva) o non abbastanza negativamente per condurre.
        if v_be >= -0.5: # Ad esempio, per un PNP al silicio si accende a circa -0.6V/-0.7V
            return 0.0

        # Corrente di base (stimata dal diodo BE)
        # La corrente Ib in un PNP scorre fuori dalla base.
        # Id_diode = self.Is * (np.exp(abs(v_be) / self.Vt) - 1.0)
        # Se v_be è negativa (es. -0.7V), abs(v_be) la rende positiva per l'esponenziale.
        # Questa corrente Ic è la corrente del diodo BE moltiplicata per Beta_F
        # e modificata dall'effetto Early e dalla regione.

        # Effetto Early sulla tensione Collettore-Emettitore
        early_factor = (1.0 - abs(v_ce) / self.Va) if self.Va > 0 else 1.0
        
        # Corrente di collettore in regione attiva inversa (per PNP è la regione attiva diretta)
        # Ic = -Beta_F * Ib (dove Ib è la corrente che entra nella base, quindi è positiva)
        # Questo è il modello semplificato, dove Ic = -Is * exp(|Vbe|/Vt) * (1 - |Vce|/Va)
        
        # Calcoliamo la corrente come valore assoluto e poi applichiamo il segno negativo
        # per indicare che esce dal collettore.
        
        # Regione Attiva: V_BE < 0, V_CE < 0, e |V_CE| > |V_BE| - V_saturation (es. 0.2V)
        # E.g. se V_BE=-0.7V, V_CE=-5V, allora |V_CE|=5 > 0.7-0.2=0.5. È in attiva.
        # Saturazione: V_BE < 0, V_CE < 0, e |V_CE| <= |V_BE| - V_saturation
        # E.g. se V_BE=-0.7V, V_CE=-0.1V, allora |V_CE|=0.1 <= 0.7-0.2=0.5. È in saturazione (triodo per MOSFET)

        # Usiamo il modello classico per un PNP, dove Ic è negativa.
        # La formula è Is * (exp(V_EB/Vt) - 1) * (1 - V_CB/Va)
        # V_EB = V_E - V_B = -V_BE = abs(v_be)
        # V_CB = V_C - V_B
        
        # Per il modello semplificato che abbiamo usato finora (corrente di collettore dominante):
        # Corrente di collettore in regione attiva forward (PNP)
        # Ic = -Is * exp(abs(v_be)/Vt) * (1 - abs(v_ce)/Va)
        
        # Consideriamo v_be è già < 0.
        # Utilizziamo la convenzione che Ic è la corrente che scorre DENTRO il collettore.
        # Per un PNP, scorre FUORI, quindi il valore sarà negativo.
        
        # Regione Attiva: V_BE < V_onset_PNP (e.g. -0.5V), e V_CE < (V_BE + V_offset_PNP) (con offset piccolo negativo, e.g. -0.2V)
        # Qui V_BE è negativa.
        # V_CE è negativa.

        # V_EB = -v_be # Tensione Emettitore-Base (positiva)
        # V_EC = -v_ce # Tensione Emettitore-Collettore (positiva)
        
        # if V_EB < 0.5: # Cutoff
        #    return 0.0
        
        # Corrente di Drain (Id) è calcolata in base a Vgs e Vds
        # Per BJT, Ic in base a Vbe e Vce.
        
        # Questo è il modello per la regione attiva
        current_active = -self.Is * (np.exp(abs(v_be) / self.Vt) - 1.0) * early_factor
        
        # Se siamo in saturazione (equivalente al triodo per MOSFET), la formula cambia.
        # La saturazione per un PNP è quando |V_CE| è piccola, es. |V_CE| < 0.2V
        # O se V_CE > V_BE - 0.2 (per NPN) -> V_CE < V_BE + 0.2 (per PNP, con Vbe negativa)
        
        # Condizione di saturazione per PNP: V_CE è molto vicino a V_BE (cioè, |V_CE| è piccola)
        # Ovvero, |V_CE| < |V_BE| - V_C_E_SAT (es. V_C_E_SAT = 0.2V)
        # Se v_ce (negativa) è "più positiva" di v_be + 0.2V (es. v_ce = -0.1, v_be = -0.7 => -0.1 > -0.5)
        # --> è in saturazione.
        
        Vce_sat_threshold = v_be + 0.2 # Esempio di soglia di saturazione per PNP
        
        if v_ce >= Vce_sat_threshold: # Se Vce è più vicina a 0 che Vbe + 0.2 (quindi in saturazione)
            # Regione di saturazione (o "triodo" per BJT)
            # La corrente in saturazione è spesso limitata da resistenze esterne.
            # Un modello molto semplificato potrebbe essere Ic_sat = -Beta_F * Ib_sat
            # O un modello più complesso che tiene conto di Vce_sat.
            # Per ora, in un modello semplificato, potremmo approssimare con la corrente attiva
            # e lasciare che il solutore gestisca i limiti fisici con le resistenze esterne.
            # Per una maggiore precisione, si implementerebbe un modello di saturazione.
            return max(current_active, -self.Beta_F * self.Is * (np.exp(abs(v_be)/self.Vt) - 1.0
