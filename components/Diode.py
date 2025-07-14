# components/diode.py
import numpy as np
from components.component import Component

class Diode(Component):
    def __init__(self, name: str, anode_node: str, cathode_node: str, Is: float = 1e-14, N: float = 1.0, Vt: float = 0.0258):
        """
        Inizializza un diodo usando il modello di Shockley.
        Args:
            name (str): Nome univoco dell'istanza (es. "D1").
            anode_node (str): Nome del nodo dell'anodo.
            cathode_node (str): Nome del nodo del catodo.
            Is (float): Corrente di saturazione inversa.
            N (float): Fattore di idealità.
            Vt (float): Tensione termica (kT/q).
        """
        super().__init__(name, anode_node, cathode_node)
        self.Is = Is
        self.N = N
        self.Vt = Vt

    def calculate_current(self, Vd: float) -> float:
        """
        Calcola la corrente attraverso il diodo data la tensione ai suoi capi (Vd = V_anode - V_cathode).
        Usa l'equazione di Shockley.
        """
        # Evita overflow per Vd molto grandi
        if Vd / (self.N * self.Vt) > 700: # Limite approssimato per exp()
            return self.Is * (np.exp(700) - 1)
        
        return self.Is * (np.exp(Vd / (self.N * self.Vt)) - 1)

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

