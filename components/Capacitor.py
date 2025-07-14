# components/capacitor.py
import numpy as np
from components.component import Component

class Capacitor(Component):
    def __init__(self, name: str, node1: str, node2: str, capacitance: float):
        """
        Inizializza un condensatore.
        Args:
            name (str): Nome univoco dell'istanza (es. "C1").
            node1 (str): Nome del primo nodo di connessione.
            node2 (str): Nome del secondo nodo di connessione.
            capacitance (float): Valore della capacità in Farad.
        """
        super().__init__(name, node1, node2)
        self.capacitance = capacitance
        # Variabili di stato per l'integrazione trapezoidale
        self.v_prev = 0.0 # Tensione ai capi del condensatore al passo precedente
        self.i_prev = 0.0 # Corrente attraverso il condensatore al passo precedente

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi del condensatore alla matrice MNA (stamp_A) e al vettore RHS (stamp_B)
        usando il metodo trapezoidale.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        node1_id, node2_id = self.node_ids

        # Conduttanza equivalente per il metodo trapezoidale
        G_eq = 2.0 * self.capacitance / dt

        # Contributi alla matrice MNA (parte dipendente dalla tensione attuale)
        if node1_id != 0: stamp_A[node1_id, node1_id] += G_eq
        if node2_id != 0: stamp_A[node2_id, node2_id] += G_eq
        if node1_id != 0 and node2_id != 0:
            stamp_A[node1_id, node2_id] -= G_eq
            stamp_A[node2_id, node1_id] -= G_eq

        # Contributi al vettore RHS (parte dipendente dallo stato precedente)
        # I_eq = G_eq * V_C_prev + I_C_prev
        # Questa corrente viene sottratta al nodo positivo e aggiunta al nodo negativo
        V_C_prev = prev_solution[node1_id] - prev_solution[node2_id] # Tensione ai capi al passo precedente
        i_eq = G_eq * V_C_prev + self.i_prev # Corrente equivalente

        if node1_id != 0: stamp_B[node1_id] -= i_eq
        if node2_id != 0: stamp_B[node2_id] += i_eq

        return stamp_A, stamp_B

    def update_state(self, v_curr: float, i_curr: float):
        """
        Aggiorna lo stato interno del condensatore per il prossimo passo temporale.
        Args:
            v_curr (float): Tensione ai capi del condensatore al passo attuale.
            i_curr (float): Corrente attraverso il condensatore al passo attuale.
        """
        self.v_prev = v_curr
        self.i_prev = i_curr

