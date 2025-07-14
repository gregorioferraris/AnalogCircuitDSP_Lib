# components/resistor.py
import numpy as np
from components.component import Component

class Resistor(Component):
    def __init__(self, name: str, node1: str, node2: str, resistance: float):
        """
        Inizializza un resistore.
        Args:
            name (str): Nome univoco dell'istanza (es. "R1").
            node1 (str): Nome del primo nodo di connessione.
            node2 (str): Nome del secondo nodo di connessione.
            resistance (float): Valore della resistenza in Ohm.
        """
        super().__init__(name, node1, node2)
        self.resistance = resistance

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi del resistore alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        node1_id, node2_id = self.node_ids # Accedi agli ID numerici dei nodi

        G = 1.0 / self.resistance

        # Contributi alla matrice MNA (ammettenze)
        # Assumiamo che il nodo 0 sia il ground e non lo gestiamo esplicitamente qui,
        # ma il solutore principale fisserà la sua tensione a 0.
        if node1_id != 0: stamp_A[node1_id, node1_id] += G
        if node2_id != 0: stamp_A[node2_id, node2_id] += G
        if node1_id != 0 and node2_id != 0:
            stamp_A[node1_id, node2_id] -= G
            stamp_A[node2_id, node1_id] -= G

        return stamp_A, stamp_B

