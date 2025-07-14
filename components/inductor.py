# components/inductor.py
import numpy as np
from components.component import Component

class Inductor(Component):
    def __init__(self, name: str, node1: str, node2: str, inductance: float):
        """
        Inizializza un induttore.
        Args:
            name (str): Nome univoco dell'istanza (es. "L1").
            node1 (str): Nome del primo nodo di connessione.
            node2 (str): Nome del secondo nodo di connessione.
            inductance (float): Valore dell'induttanza in Henry.
        """
        super().__init__(name, node1, node2)
        self.L = inductance
        # Variabili di stato per l'integrazione trapezoidale
        self.i_prev = 0.0 # Corrente attraverso l'induttore al passo precedente
        self.v_prev = 0.0 # Tensione ai capi dell'induttore al passo precedente

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi dell'induttore alla matrice MNA (stamp_A) e al vettore RHS (stamp_B)
        usando il metodo trapezoidale.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        node1_id, node2_id = self.node_ids

        # Resistenza equivalente per il metodo trapezoidale
        # R_eq = 2L / dt
        # G_eq = dt / (2L)
        G_eq = dt / (2.0 * self.L)

        # Contributi alla matrice MNA (ammettenze)
        if node1_id != 0: stamp_A[node1_id, node1_id] += G_eq
        if node2_id != 0: stamp_A[node2_id, node2_id] += G_eq
        if node1_id != 0 and node2_id != 0:
            stamp_A[node1_id, node2_id] -= G_eq
            stamp_A[node2_id, node1_id] -= G_eq

        # Contributi al vettore RHS (parte dipendente dallo stato precedente)
        # V_eq = I_L_prev * (2L / dt) + V_L_prev
        # Questa tensione viene applicata come sorgente di tensione equivalente
        V_eq = self.i_prev * (2.0 * self.L / dt) + self.v_prev

        # L'induttore è come una sorgente di tensione in serie con una resistenza
        # La corrente che esce da node1 e entra in node2 è (V_node1 - V_node2 - V_eq) / (2L/dt)
        # O, in termini di MNA, aggiungiamo una corrente equivalente al vettore B
        # I_eq = G_eq * V_eq
        i_eq = G_eq * V_eq

        if node1_id != 0: stamp_B[node1_id] -= i_eq
        if node2_id != 0: stamp_B[node2_id] += i_eq

        return stamp_A, stamp_B

    def update_state(self, v_curr: float, i_curr: float):
        """
        Aggiorna lo stato interno dell'induttore per il prossimo passo temporale.
        Args:
            v_curr (float): Tensione ai capi dell'induttore al passo attuale.
            i_curr (float): Corrente attraverso l'induttore al passo attuale.
        """
        self.v_prev = v_curr
        self.i_prev = i_curr

