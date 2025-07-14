# components/schottky_diode.py
import numpy as np
from components.component import Component

class SchottkyDiode(Component):
    def __init__(self, name: str, anode_node: str, cathode_node: str, Is: float = 1e-9, N: float = 1.05, Vt: float = 0.0258):
        """
        Inizializza un diodo Schottky.
        Simile a un diodo standard, ma con tensione di forward più bassa e recupero inverso più veloce.
        Args:
            name (str): Nome univoco dell'istanza (es. "SD1").
            anode_node (str): Nome del nodo dell'anodo.
            cathode_node (str): Nome del nodo del catodo.
            Is (float): Corrente di saturazione inversa (tipicamente più alta dei diodi PN).
            N (float): Fattore di idealità (spesso vicino a 1).
            Vt (float): Tensione termica (kT/q).
        """
        super().__init__(name, anode_node, cathode_node)
        self.Is = Is
        self.N = N
        self.Vt = Vt

    def calculate_current(self, Vd: float) -> float:
        """
        Calcola la corrente attraverso il diodo Schottky data la tensione ai suoi capi (Vd = V_anode - V_cathode).
        Usa l'equazione di Shockley.
        """
        if Vd / (self.N * self.Vt) > 700: # Limite per evitare overflow
            return self.Is * (np.exp(700) - 1)
        
        current = self.Is * (np.exp(Vd / (self.N * self.Vt)) - 1)
        
        # Assicurati che la corrente non sia negativa in polarizzazione inversa (idealmente)
        return max(0.0, current)

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

