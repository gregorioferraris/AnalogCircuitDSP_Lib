# components/jfet.py
import numpy as np
from components.component import Component

class JFET(Component):
    def __init__(self, name: str, drain_node: str, gate_node: str, source_node: str,
                 Idss: float = 0.01, Vp: float = -2.0, type: str = 'nfet'):
        """
        Inizializza un JFET (modello semplificato).
        Args:
            name (str): Nome univoco dell'istanza (es. "J1").
            drain_node (str): Nome del nodo di drain.
            gate_node (str): Nome del nodo di gate.
            source_node (str): Nome del nodo di source.
            Idss (float): Corrente di drain a Vgs=0 e Vds > Vp (corrente di saturazione).
            Vp (float): Tensione di pinch-off (per NFET è negativa, per PFET è positiva).
            type (str): Tipo di JFET ('nfet' o 'pfet').
        """
        super().__init__(name, drain_node, gate_node, source_node)
        self.pin_names = ('drain', 'gate', 'source') # Per mappatura nodi con nome
        self.Idss = Idss
        self.Vp = Vp
        self.type = type.lower()
        if self.type not in ['nfet', 'pfet']:
            raise ValueError("JFET type must be 'nfet' or 'pfet'.")

    def calculate_drain_current(self, Vgs: float, Vds: float) -> float:
        """
        Calcola la corrente di drain (Id) per un JFET (modello base).
        Args:
            Vgs (float): Tensione Gate-Source.
            Vds (float): Tensione Drain-Source.
        Returns:
            float: Corrente di drain.
        """
        Id = 0.0
        if self.type == 'nfet':
            if Vgs >= self.Vp: # Condizione per la conduzione (Vp è negativo per NFET)
                if Vds > (Vgs - self.Vp): # Saturazione
                    Id = self.Idss * (1 - Vgs / self.Vp)**2
                else: # Regione ohmica (triode)
                    # Modello semplificato per la regione ohmica
                    Id = self.Idss * (1 - Vgs / self.Vp)**2 * (1 - (Vds / (Vgs - self.Vp) - 1)**2)
                    # Un modello più comune per la regione ohmica è:
                    # Id = self.Idss * (2 * (Vgs - self.Vp) * Vds - Vds**2) / self.Vp**2
                    # Ma richiede un'implementazione più attenta per la transizione.
                    # Per ora, usiamo una forma che assicura Id=0 se Vds=0 e Vgs=Vp
                    Id = self.Idss * (1 - Vgs / self.Vp)**2 * (Vds / (Vgs - self.Vp)) # Semplificazione
                    Id = np.clip(Id, 0, self.Idss) # Limita la corrente
            else: # Cutoff
                Id = 0.0
        elif self.type == 'pfet':
            # Per PFET, le tensioni sono invertite (Vsg, Vsd) e Vp è positivo
            Vsg = -Vgs
            Vsd = -Vds
            if Vsg >= -self.Vp: # Condizione per la conduzione (Vp è positivo per PFET)
                if Vsd > (Vsg - (-self.Vp)): # Saturazione
                    Id = self.Idss * (1 - Vsg / (-self.Vp))**2
                else: # Regione ohmica (triode)
                    Id = self.Idss * (1 - Vsg / (-self.Vp))**2 * (Vsd / (Vsg - (-self.Vp)))
                    Id = np.clip(Id, 0, self.Idss)
            else: # Cutoff
                Id = 0.0
            Id = -Id # Corrente di drain per PFET è in direzione opposta

        return Id

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

