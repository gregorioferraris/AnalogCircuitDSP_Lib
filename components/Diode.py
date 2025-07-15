import numpy as np
from components.component import Component

class Diode(Component):
    def __init__(self, name: str, anode_node: str, cathode_node: str, Is: float = 1e-14, N: float = 1.0, Vt: float = 0.0258):
        """
        Inizializza un diodo a giunzione PN.
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
        # Aggiungo un piccolo limite per evitare overflow nell'esponenziale
        arg = Vd / (self.N * self.Vt)
        if arg > 700: # Per prevenire overflow numerico
            return self.Is * (np.exp(700) - 1)
        elif arg < -70: # Per prevenire underflow e risultati inaccurati in inversa molto profonda
            return -self.Is
        else:
            return self.Is * (np.exp(arg) - 1)

    def calculate_conductance(self, Vd: float) -> float:
        """
        Calcola la conduttanza dinamica del diodo (dI/dVd).
        Necessaria per il solutore Newton-Raphson.
        """
        arg = Vd / (self.N * self.Vt)
        if arg > 700: # Per la stabilità numerica in avanti
            # La conduttanza diventa molto alta in conduzione forte, usiamo un limite
            return self.Is / (self.N * self.Vt) * np.exp(700)
        elif arg < -70: # La conduttanza è molto piccola in inversa
            return 1e-9 # Piccola conduttanza per stabilità numerica, non 0
        else:
            return self.Is / (self.N * self.Vt) * np.exp(arg)

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

