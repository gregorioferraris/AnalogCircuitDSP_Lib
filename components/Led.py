# components/led.py
import numpy as np
from components.component import Component

class LED(Component):
    def __init__(self, name: str, anode_node: str, cathode_node: str,
                 Is: float = 1e-12, N: float = 1.8, Vt: float = 0.0258, V_fwd_threshold: float = 1.8):
        """
        Inizializza un LED (Light Emitting Diode) usando un modello di diodo con soglia.
        Args:
            name (str): Nome univoco dell'istanza (es. "LED1").
            anode_node (str): Nome del nodo dell'anodo.
            cathode_node (str): Nome del nodo del catodo.
            Is (float): Corrente di saturazione inversa.
            N (float): Fattore di idealità.
            Vt (float): Tensione termica (kT/q).
            V_fwd_threshold (float): Tensione di soglia di conduzione (es. 1.8V per rosso, 3V per blu).
        """
        super().__init__(name, anode_node, cathode_node)
        self.Is = Is
        self.N = N
        self.Vt = Vt
        self.V_fwd_threshold = V_fwd_threshold # Tensione di soglia specifica del LED

    def calculate_current(self, Vd: float) -> float:
        """
        Calcola la corrente attraverso il LED data la tensione ai suoi capi (Vd = V_anode - V_cathode).
        Modello di Shockley con una soglia di tensione più pronunciata.
        """
        # Per modellare la soglia, possiamo sottrarre la V_fwd_threshold dalla tensione applicata
        # o usare una funzione che simuli una conduzione molto bassa sotto soglia.
        # Qui usiamo un approccio che sposta la caratteristica Shockley.
        
        Vd_effective = Vd - self.V_fwd_threshold
        
        if Vd_effective / (self.N * self.Vt) > 700: # Limite per evitare overflow
            return self.Is * (np.exp(700) - 1)
        
        current = self.Is * (np.exp(Vd_effective / (self.N * self.Vt)) - 1)
        
        # Assicurati che la corrente non sia negativa in polarizzazione inversa
        return max(0.0, current)

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

