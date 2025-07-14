# components/mosfet.py
import numpy as np
from components.component import Component

class MOSFET(Component):
    def __init__(self, name: str, drain_node: str, gate_node: str, source_node: str,
                 Vth: float = 0.5, Kp: float = 1e-3, type: str = 'nmos'):
        """
        Inizializza un MOSFET (modello semplificato).
        Args:
            name (str): Nome univoco dell'istanza (es. "M1").
            drain_node (str): Nome del nodo di drain.
            gate_node (str): Nome del nodo di gate.
            source_node (str): Nome del nodo di source.
            Vth (float): Tensione di soglia.
            Kp (float): Fattore di transconduttanza (mu_n * C_ox * W/L).
            type (str): Tipo di MOSFET ('nmos' o 'pmos').
        """
        super().__init__(name, drain_node, gate_node, source_node)
        self.pin_names = ('drain', 'gate', 'source') # Per mappatura nodi con nome
        self.Vth = Vth
        self.Kp = Kp
        self.type = type.lower()
        if self.type not in ['nmos', 'pmos']:
            raise ValueError("MOSFET type must be 'nmos' or 'pmos'.")

    def calculate_drain_current(self, Vgs: float, Vds: float) -> float:
        """
        Calcola la corrente di drain (Id) per un MOSFET (modello base).
        Args:
            Vgs (float): Tensione Gate-Source.
            Vds (float): Tensione Drain-Source.
        Returns:
            float: Corrente di drain.
        """
        Id = 0.0
        if self.type == 'nmos':
            if Vgs < self.Vth: # Cutoff
                Id = 0.0
            elif Vds < (Vgs - self.Vth): # Triode (Lineare)
                Id = self.Kp * ((Vgs - self.Vth) * Vds - 0.5 * Vds**2)
            else: # Saturazione
                Id = 0.5 * self.Kp * (Vgs - self.Vth)**2
        elif self.type == 'pmos':
            # Per PMOS, le tensioni sono invertite (Vsg, Vsd)
            Vsg = -Vgs
            Vsd = -Vds
            if Vsg < -self.Vth: # Cutoff (Vth per PMOS è negativo)
                Id = 0.0
            elif Vsd < (Vsg - (-self.Vth)): # Triode (Lineare)
                Id = self.Kp * ((Vsg - (-self.Vth)) * Vsd - 0.5 * Vsd**2)
            else: # Saturazione
                Id = 0.5 * self.Kp * (Vsg - (-self.Vth))**2
            Id = -Id # Corrente di drain per PMOS è in direzione opposta

        return Id

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

