import numpy as np
from components.component import Component

class Splitter(Component):
    def __init__(self, name: str, input_node: str, *output_nodes: str):
        """
        Inizializza un componente Splitter di segnale.
        Prende un segnale da un nodo di ingresso e lo replica su più nodi di uscita.
        Questo modello impone che le tensioni di tutti i nodi di uscita siano uguali
        alla tensione del nodo di ingresso.
        Args:
            name (str): Nome univoco dell'istanza (es. "SPLITTER1").
            input_node (str): Nome del nodo di ingresso.
            *output_nodes (str): Nomi dei nodi di uscita (uno o più).
        """
        # I nodi del componente saranno (input_node, output_node_1, output_node_2, ...)
        super().__init__(name, input_node, *output_nodes)
        self.pin_names = ('input',) + tuple(f'output_{i+1}' for i in range(len(output_nodes)))
        
        self.input_node_name = input_node
        self.output_node_names = output_nodes
        self.num_outputs = len(output_nodes)

        # Indici per le correnti ausiliarie aggiunte al sistema MNA
        # Ogni output_node avrà una sua corrente ausiliaria per imporre la relazione di tensione
        self.output_current_indices = [-1] * self.num_outputs

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi dello Splitter alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        Ogni nodo di uscita aggiunge una riga e una colonna extra al sistema MNA
        per imporre la relazione V_output = V_input.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        # Recupera gli ID dei nodi
        # node_ids è una tupla: (input_id, output_1_id, output_2_id, ...)
        input_id = self.node_ids[0]
        output_ids = self.node_ids[1:] # Tutti gli ID dei nodi di uscita

        # Per ogni nodo di uscita, aggiungiamo un'equazione ausiliaria
        for i, output_id in enumerate(output_ids):
            aux_current_idx = self.output_current_indices[i] # Indice della corrente ausiliaria per questo output

            # --- Equazioni KCL (correnti ausiliarie) ---
            # La corrente ausiliaria entra nel nodo di uscita ed esce dal nodo di ingresso
            # (o viceversa, la convenzione qui è che impone V_out - V_in = 0)
            
            # Contributo della corrente ausiliaria ai nodi di uscita e ingresso
            # Questa corrente è la "corrente" che impone il cortocircuito concettuale
            if output_id != 0: stamp_A[output_id, aux_current_idx] += 1
            if input_id != 0: stamp_A[input_id, aux_current_idx] -= 1 # La corrente "esce" dall'input per imporre l'uguaglianza

            # --- Equazione ausiliaria: V_output - V_input = 0 ---
            # Questa equazione viene aggiunta alla riga MNA corrispondente a aux_current_idx
            if output_id != 0: stamp_A[aux_current_idx, output_id] += 1
            if input_id != 0: stamp_A[aux_current_idx, input_id] -= 1
            stamp_B[aux_current_idx] = 0.0 # Il lato destro è zero per V_output = V_input

        return stamp_A, stamp_B

    def _set_output_current_indices(self, indices: list[int]):
        """Metodo interno chiamato dal solutore per assegnare gli indici delle correnti ausiliarie di uscita."""
        if len(indices) != self.num_outputs:
            raise ValueError(f"Numero di indici forniti ({len(indices)}) non corrisponde al numero di uscite ({self.num_outputs}).")
        self.output_current_indices = indices

