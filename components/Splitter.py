# components/component.py

import numpy as np

class Component:
    """
    Classe base per tutti i componenti del circuito.
    Definisce l'interfaccia comune per la gestione dei nodi e l'identificazione.
    Include la logica per l'applicazione delle tolleranze di produzione.
    """
    def __init__(self, name: str, *node_names_str: str):
        """
        Inizializza un componente.
        Args:
            name (str): Nome univoco dell'istanza del componente (es. "R1", "C_bypass").
            *node_names_str (str): Nomi dei nodi a cui il componente è connesso (stringhe, es. "input", "output").
        """
        self.name = name
        self.node_names_str = node_names_str # Nomi dei nodi come stringhe (es. ('input', 'output'))
        self.node_ids = None # Verrà popolato dalla classe Circuit con gli ID numerici dei nodi
        self.component_id = None # Verrà assegnato dalla classe Circuit

    def _set_node_ids(self, ids):
        """
        Metodo interno chiamato dalla classe Circuit per impostare gli ID numerici dei nodi.
        'ids' può essere una tupla (ID_NODO1, ID_NODO2, ...) o un dizionario {'nome_pin': ID_NODO}.
        """
        self.node_ids = ids

    def _set_component_id(self, comp_id):
        """Metodo interno chiamato dalla classe Circuit per assegnare un ID univoco al componente."""
        self.component_id = comp_id

    def _apply_tolerance(self, nominal_value: float, tolerance_percent: float, distribution: str = 'uniform') -> float:
        """
        Applica una tolleranza percentuale al valore nominale.
        Args:
            nominal_value (float): Il valore nominale del parametro.
            tolerance_percent (float): La tolleranza in percentuale (es. 5.0 per +/- 5%).
            distribution (str): Tipo di distribuzione ('uniform' o 'normal').
        Returns:
            float: Il valore del parametro con la tolleranza applicata.
        """
        if tolerance_percent == 0:
            return nominal_value
        
        if distribution == 'uniform':
            # Genera un fattore casuale tra -tolerance_percent e +tolerance_percent
            variation_factor = (2 * np.random.rand() - 1) * (tolerance_percent / 100.0)
            return nominal_value * (1 + variation_factor)
        elif distribution == 'normal':
            # Assumiamo che la tolleranza percentuale rappresenti 3 deviazioni standard (3-sigma)
            # Questo copre circa il 99.7% dei valori
            std_dev = nominal_value * (tolerance_percent / (3 * 100.0))
            return np.random.normal(nominal_value, std_dev)
        else:
            raise ValueError("Tipo di distribuzione non supportato. Usare 'uniform' o 'normal'.")

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Metodo astratto per ottenere i "contributi" del componente alla matrice MNA e al vettore RHS.
        Deve essere implementato dalle sottoclassi.
        Args:
            num_total_equations (int): Dimensione totale del sistema di equazioni.
            dt (float): Passo temporale (1/Fs).
            current_solution_guess (np.ndarray): Guess corrente per le incognite (tensioni ai nodi, correnti sorgenti V).
            prev_solution (np.ndarray): Soluzione del passo temporale precedente.
            time (float): Tempo attuale della simulazione.
        Returns:
            tuple: (stamp_A, stamp_B) - Matrice e vettore dei contributi del componente.
        """
        raise NotImplementedError("Il metodo get_stamps deve essere implementato dalle sottoclassi.")

    def __str__(self):
        return f"{self.__class__.__name__}(Name: '{self.name}', Nodes: {self.node_names_str})"

