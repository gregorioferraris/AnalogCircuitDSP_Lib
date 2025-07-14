# components/ldr.py
import numpy as np
from components.component import Component

class LDR(Component):
    def __init__(self, name: str, node1: str, node2: str, dark_resistance: float = 1e6, light_resistance: float = 1e3, light_level: float = 0.0):
        """
        Inizializza una LDR (Light Dependent Resistor).
        La sua resistenza varia in base al livello di luce.
        Args:
            name (str): Nome univoco dell'istanza (es. "LDR1").
            node1 (str): Nome del primo nodo di connessione.
            node2 (str): Nome del secondo nodo di connessione.
            dark_resistance (float): Resistenza in assenza di luce (Ohm).
            light_resistance (float): Resistenza in piena luce (Ohm).
            light_level (float): Livello di luce attuale (0.0 = buio, 1.0 = piena luce).
        """
        super().__init__(name, node1, node2)
        self.dark_resistance = dark_resistance
        self.light_resistance = light_resistance
        self._light_level = light_level # Livello di luce attuale (0.0 a 1.0)

    def set_light_level(self, level: float):
        """
        Imposta il livello di luce per la LDR (tra 0.0 e 1.0).
        """
        self._light_level = np.clip(level, 0.0, 1.0)

    def get_resistance(self) -> float:
        """
        Calcola la resistenza attuale della LDR in base al livello di luce.
        Modello semplificato: interpolazione logaritmica tra resistenza al buio e alla luce.
        """
        # Interpolazione logaritmica per una variazione più realistica
        log_dark_R = np.log10(self.dark_resistance)
        log_light_R = np.log10(self.light_resistance)
        
        log_current_R = log_dark_R + (log_light_R - log_dark_R) * self._light_level
        return 10**log_current_R

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi della LDR alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        La LDR si comporta come un resistore, ma il suo valore di resistenza è dinamico.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        node1_id, node2_id = self.node_ids

        G = 1.0 / self.get_resistance() # Usa la resistenza attuale

        if node1_id != 0: stamp_A[node1_id, node1_id] += G
        if node2_id != 0: stamp_A[node2_id, node2_id] += G
        if node1_id != 0 and node2_id != 0:
            stamp_A[node1_id, node2_id] -= G
            stamp_A[node2_id, node1_id] -= G

        return stamp_A, stamp_B

