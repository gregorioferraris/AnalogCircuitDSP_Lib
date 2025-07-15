# circuit_solver/circuit.py (AGGIORNATO)

import numpy as np
from components.component import Component
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor
from components.voltage_source import VoltageSource
from components.current_source import CurrentSource
from components.diode import Diode
from components.transformer import Transformer
from components.splitter import Splitter # Nuovo import

class Circuit:
    def __init__(self):
        self.components = []
        self.node_names = ['0'] # Nodo di terra per convenzione
        self.node_map = {'0': 0}
        self.num_nodes = 1
        self.voltage_source_count = 0 # Contatore per gli indici delle correnti delle Vs
        self.transformer_count = 0    # Contatore per gli indici delle correnti dei trasformatori
        self.splitter_current_count = 0 # Nuovo contatore per gli indici delle correnti degli splitter

    def add_component(self, component: Component):
        self.components.append(component)
        
        # Registra i nomi dei nodi e assegna ID numerici se non già presenti
        for node_name in component.node_names:
            if node_name not in self.node_map:
                self.node_map[node_name] = self.num_nodes
                self.node_names.append(node_name)
                self.num_nodes += 1
        
        # Assegna gli ID numerici ai nodi del componente
        # node_names è una tupla di stringhe, node_ids sarà una tupla di int
        component.node_ids = tuple(self.node_map[name] for name in component.node_names)

        # Pre-assegna gli indici per le correnti delle sorgenti di tensione
        if isinstance(component, VoltageSource):
            component._set_current_index(self.num_nodes + self.voltage_source_count)
            self.voltage_source_count += 1
        
        # Pre-assegna gli indici per le correnti dei trasformatori
        elif isinstance(component, Transformer):
            # Il trasformatore aggiunge DUE correnti ausiliarie
            ip_idx = self.num_nodes + self.voltage_source_count + (self.transformer_count * 2)
            is_idx = ip_idx + 1
            component._set_current_indices(ip_idx, is_idx)
            self.transformer_count += 1
        
        # Pre-assegna gli indici per le correnti dello splitter
        elif isinstance(component, Splitter):
            # Lo splitter aggiunge una corrente ausiliaria per ogni nodo di uscita
            start_idx = self.num_nodes + self.voltage_source_count + (self.transformer_count * 2) + self.splitter_current_count
            indices = [start_idx + i for i in range(component.num_outputs)]
            component._set_output_current_indices(indices)
            self.splitter_current_count += component.num_outputs


    def get_num_total_equations(self):
        """Restituisce il numero totale di equazioni (nodi + correnti Vs + correnti Trasformatore + correnti Splitter)."""
        return self.num_nodes + self.voltage_source_count + (self.transformer_count * 2) + self.splitter_current_count

