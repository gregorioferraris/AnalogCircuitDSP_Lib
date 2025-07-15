import numpy as np
from components.component import Component

class Transformer(Component):
    def __init__(self, name: str, prim_plus_node: str, prim_minus_node: str,
                 sec_plus_node: str, sec_minus_node: str,
                 turns_ratio: float = 1.0):
        """
        Inizializza un modello di trasformatore ideale.
        Args:
            name (str): Nome univoco dell'istanza (es. "XFMR1").
            prim_plus_node (str): Nome del nodo positivo del primario.
            prim_minus_node (str): Nome del nodo negativo del primario.
            sec_plus_node (str): Nome del nodo positivo del secondario.
            sec_minus_node (str): Nome del nodo negativo del secondario.
            turns_ratio (float): Rapporto di spire N_sec / N_prim.
                                 Se turns_ratio > 1, è step-up.
                                 Se turns_ratio < 1, è step-down.
        """
        super().__init__(name, prim_plus_node, prim_minus_node, sec_plus_node, sec_minus_node)
        self.pin_names = ('prim_plus', 'prim_minus', 'sec_plus', 'sec_minus')

        self.turns_ratio = turns_ratio # N_sec / N_prim

        # Indici per le correnti primarie e secondarie aggiunte al sistema MNA
        self.ip_index = -1 # Indice della corrente che attraversa il primario
        self.is_index = -1 # Indice della corrente che attraversa il secondario

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi del trasformatore ideale alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        Per un trasformatore ideale, si aggiungono due equazioni e due incognite di corrente.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        # Recupera gli ID dei nodi
        # node_ids è una tupla: (prim_plus_id, prim_minus_id, sec_plus_id, sec_minus_id)
        prim_plus_id = self.node_ids[0]
        prim_minus_id = self.node_ids[1]
        sec_plus_id = self.node_ids[2]
        sec_minus_id = self.node_ids[3]
        
        # Recupera gli indici delle correnti del trasformatore
        ip_idx = self.ip_index # Corrente nel primario (variabile ausiliaria)
        is_idx = self.is_index # Corrente nel secondario (variabile ausiliaria)

        # --- Equazioni KCL (correnti che attraversano il primario e il secondario) ---
        # Si assume che ip_idx sia la corrente che entra nel prim_plus ed esce dal prim_minus
        # E is_idx sia la corrente che entra nel sec_plus ed esce dal sec_minus

        # Contributo della corrente del primario ai nodi prim_plus e prim_minus
        if prim_plus_id != 0: stamp_A[prim_plus_id, ip_idx] += 1 # La corrente ip entra in prim_plus
        if prim_minus_id != 0: stamp_A[prim_minus_id, ip_idx] -= 1 # La corrente ip esce da prim_minus

        # Contributo della corrente del secondario ai nodi sec_plus e sec_minus
        if sec_plus_id != 0: stamp_A[sec_plus_id, is_idx] += 1 # La corrente is entra in sec_plus
        if sec_minus_id != 0: stamp_A[sec_minus_id, is_idx] -= 1 # La corrente is esce da sec_minus

        # --- Equazioni aggiuntive (regole del trasformatore ideale) ---
        # 1. Relazione tensione: V_sec = V_prim * (N_sec / N_prim)
        #    (V_sec_plus - V_sec_minus) - turns_ratio * (V_prim_plus - V_prim_minus) = 0
        # Questa equazione è posta nella riga di MNA associata a ip_idx
        stamp_A[ip_idx, sec_plus_id] += 1
        stamp_A[ip_idx, sec_minus_id] -= 1
        stamp_A[ip_idx, prim_plus_id] -= self.turns_ratio
        stamp_A[ip_idx, prim_minus_id] += self.turns_ratio
        stamp_B[ip_idx] = 0.0 # Il lato destro è zero per un trasformatore ideale

        # 2. Relazione corrente: I_prim * N_prim = I_sec * N_sec => I_prim = I_sec * (N_sec / N_prim)
        #    I_prim - turns_ratio * I_sec = 0
        # Questa equazione è posta nella riga di MNA associata a is_idx
        stamp_A[is_idx, ip_idx] += 1
        stamp_A[is_idx, is_idx] -= self.turns_ratio
        stamp_B[is_idx] = 0.0

        return stamp_A, stamp_B

    def _set_current_indices(self, ip_index: int, is_index: int):
        """Metodo interno chiamato dal solutore per assegnare gli indici delle correnti del trasformatore."""
        self.ip_index = ip_index
        self.is_index = is_index
