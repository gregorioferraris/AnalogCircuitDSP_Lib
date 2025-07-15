import numpy as np
from components.component import Component

class CurrentSource(Component):
    def __init__(self, name: str, node_plus: str, node_minus: str, initial_current: float = 0.0):
        """
        Inizializza una sorgente di corrente indipendente.
        Args:
            name (str): Nome univoco dell'istanza (es. "I_bias").
            node_plus (str): Nodo da cui la corrente esce (convenzionale).
            node_minus (str): Nodo in cui la corrente entra.
            initial_current (float): Valore iniziale della corrente (per sorgenti DC o come base).
        """
        super().__init__(name, node_plus, node_minus)
        self._current_value = initial_current

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi della sorgente di corrente alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        Una sorgente di corrente aggiunge termini solo al vettore RHS.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        node_plus_id, node_minus_id = self.node_ids
        
        current_val = self.get_current(time)

        # La corrente esce da node_plus_id ed entra in node_minus_id
        if node_plus_id != 0: stamp_B[node_plus_id] -= current_val
        if node_minus_id != 0: stamp_B[node_minus_id] += current_val

        return stamp_A, stamp_B

    def get_current(self, time: float) -> float:
        """
        Restituisce il valore della corrente della sorgente al tempo specificato.
        Questo metodo può essere sovrascritto per definire segnali dipendenti dal tempo (AC, impulso, ecc.).
        """
        return self._current_value

    def set_current(self, value: float):
        """
        Imposta il valore della corrente della sorgente.
        """
        self._current_value = value
