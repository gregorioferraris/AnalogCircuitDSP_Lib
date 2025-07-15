import numpy as np
from circuit_solver.circuit import Circuit
from components.component import Component

class Subcircuit(Circuit):
    """
    Classe base per un blocco di circuito riutilizzabile.
    Permette di definire esplicitamente i nodi di ingresso e uscita del blocco.
    """
    def __init__(self, name: str, input_nodes: list[str] = None, output_nodes: list[str] = None):
        """
        Inizializza un sottocircuito.
        Args:
            name (str): Nome univoco dell'istanza del sottocircuito.
            input_nodes (list[str]): Nomi dei nodi che fungono da ingresso per questo blocco.
            output_nodes (list[str]): Nomi dei nodi che fungono da uscita per questo blocco.
        """
        super().__init__() # Inizializza la parte Circuit
        self.name = name
        self._input_nodes_names = input_nodes if input_nodes is not None else []
        self._output_nodes_names = output_nodes if output_nodes is not None else []

        # Assicurati che i nodi di I/O siano aggiunti al circuito interno del sottocircuito
        for node_name in self._input_nodes_names:
            self.add_node(node_name)
        for node_name in self._output_nodes_names:
            self.add_node(node_name)

    def get_input_node_ids(self) -> list[int]:
        """Restituisce gli ID numerici dei nodi di ingresso del sottocircuito."""
        return [self.node_map[name] for name in self._input_nodes_names]

    def get_output_node_ids(self) -> list[int]:
        """Restituisce gli ID numerici dei nodi di uscita del sottocircuito."""
        return [self.node_map[name] for name in self._output_nodes_names]

    def get_input_node_names(self) -> list[str]:
        """Restituisce i nomi delle stringhe dei nodi di ingresso del sottocircuito."""
        return self._input_nodes_names

    def get_output_node_names(self) -> list[str]:
        """Restituisce i nomi delle stringhe dei nodi di uscita del sottocircuito."""
        return self._output_nodes_names

    def __str__(self):
        return (f"Subcircuit(Name: '{self.name}', "
                f"Inputs: {self._input_nodes_names}, "
                f"Outputs: {self._output_nodes_names}, "
                f"Components: {len(self.components)})")

