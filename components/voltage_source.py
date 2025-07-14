# components/voltage_source.py
import numpy as np
from components.component import Component

class VoltageSource(Component):
    def __init__(self, name: str, node_plus: str, node_minus: str, initial_voltage: float = 0.0):
        """
        Inizializza una sorgente di tensione indipendente.
        Args:
            name (str): Nome univoco dell'istanza (es. "V_in", "VCC").
            node_plus (str): Nome del nodo positivo.
            node_minus (str): Nome del nodo negativo.
            initial_voltage (float): Valore iniziale della tensione (per sorgenti DC o come base).
        """
        super().__init__(name, node_plus, node_minus)
        self._voltage_value = initial_voltage # Valore attuale della sorgente
        self.current_index = -1 # Indice della riga/colonna per la corrente della sorgente nel sistema MNA

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi della sorgente di tensione alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        Una sorgente di tensione aggiunge una riga e una colonna extra al sistema MNA.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        node_plus_id, node_minus_id = self.node_ids
        vs_current_idx = self.current_index # Indice della corrente della sorgente nel vettore delle incognite

        # Contributi alla KCL ai nodi (corrente della Vs)
        # La corrente esce dal nodo positivo ed entra nel nodo negativo
        if node_plus_id != 0: stamp_A[node_plus_id, vs_current_idx] += 1
        if node_minus_id != 0: stamp_A[node_minus_id, vs_current_idx] -= 1

        # Equazione della sorgente di tensione: V_node_plus - V_node_minus = V_source
        if node_plus_id != 0: stamp_A[vs_current_idx, node_plus_id] += 1
        if node_minus_id != 0: stamp_A[vs_current_idx, node_minus_id] -= 1
        
        stamp_B[vs_current_idx] = self.get_voltage(time) # Il valore della sorgente

        return stamp_A, stamp_B

    def get_voltage(self, time: float) -> float:
        """
        Restituisce il valore della tensione della sorgente al tempo specificato.
        Questo metodo può essere sovrascritto per definire segnali dipendenti dal tempo (AC, impulso, ecc.).
        """
        # Per default, restituisce il valore DC o l'ultimo valore impostato
        return self._voltage_value

    def set_voltage(self, value: float):
        """
        Imposta il valore della tensione della sorgente. Utile per sorgenti controllate esternamente.
        """
        self._voltage_value = value

    def _set_current_index(self, index: int):
        """Metodo interno chiamato dal solutore per assegnare l'indice della corrente della sorgente."""
        self.current_index = index

