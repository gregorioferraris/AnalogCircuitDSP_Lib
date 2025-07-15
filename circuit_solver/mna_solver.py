# circuit_solver/circuit.py

import numpy as np
from collections import defaultdict
from components.component import Component
from components.voltage_source import VoltageSource
from components.splitter import Splitter # Importa il nuovo Splitter

class Circuit:
    """
    Rappresenta il circuito elettrico, gestendo i nodi, i componenti
    e la mappatura tra nomi dei nodi/componenti e ID numerici.
    """
    def __init__(self):
        self.nodes = {'0': 0}  # Dizionario per mappare i nomi dei nodi ai loro ID numerici. '0' è il ground.
        self.node_map = {'0': 0} # Mappa inversa: nome -> ID
        self.next_node_id = 1  # Il prossimo ID numerico disponibile per un nodo
        self.components = []   # Lista di tutti i componenti nel circuito
        self.voltage_sources = [] # Lista delle sorgenti di tensione (aggiungono variabili ausiliarie)
        self.splitters = [] # Lista degli splitter (possono aggiungere variabili ausiliarie)
        self.subcircuits = [] # Lista dei sottocircuiti

    def add_node(self, node_name: str):
        """
        Aggiunge un nodo al circuito se non esiste già.
        Args:
            node_name (str): Il nome del nodo da aggiungere.
        Returns:
            int: L'ID numerico del nodo.
        """
        if node_name not in self.nodes:
            self.nodes[node_name] = self.next_node_id
            self.node_map[node_name] = self.next_node_id
            self.next_node_id += 1
        return self.nodes[node_name]

    def add_component(self, component: Component):
        """
        Aggiunge un componente al circuito e assegna gli ID numerici ai suoi nodi.
        Args:
            component (Component): L'istanza del componente da aggiungere.
        """
        # Assegna un ID univoco al componente
        component.component_id = len(self.components)

        # Assegna gli ID numerici ai nodi del componente
        node_ids = []
        for node_name in component.node_names_str:
            node_ids.append(self.add_node(node_name))
        component._set_node_ids(tuple(node_ids)) # Passa la tupla degli ID

        self.components.append(component)

        # Categorizza i componenti che aggiungono variabili ausiliarie
        if isinstance(component, VoltageSource):
            self.voltage_sources.append(component)
        if isinstance(component, Splitter): # Aggiungi il nuovo tipo di componente
            self.splitters.append(component)

        print(f"Aggiunto componente: {component.name} ({component.__class__.__name__})")

    def add_block(self, subcircuit: 'Subcircuit'):
        """
        Aggiunge un sottocircuito (blocco) al circuito principale.
        Args:
            subcircuit (Subcircuit): L'istanza del sottocircuito da aggiungere.
        """
        self.subcircuits.append(subcircuit)
        # Aggiungi tutti i componenti interni del sottocircuito al circuito principale
        for comp in subcircuit.components:
            self.add_component(comp)
        
        # Mappa i nodi di input/output del sottocircuito ai nodi del circuito principale
        # Questo è gestito dalla logica di Subcircuit.connect_blocks() o aggiungendo i nodi
        # direttamente al circuito principale.
        for node_name in subcircuit.nodes:
            self.add_node(node_name)

        print(f"Aggiunto sottocircuito: {subcircuit.name}")

    def get_num_nodes(self) -> int:
        """Restituisce il numero totale di nodi nel circuito (inclusi i nodi interni dei sottocircuiti)."""
        return self.next_node_id

    def get_components(self) -> list[Component]:
        """Restituisce la lista di tutti i componenti nel circuito."""
        return self.components

    def get_voltage_sources(self) -> list[VoltageSource]:
        """Restituisce la lista delle sorgenti di tensione nel circuito."""
        return self.voltage_sources
    
    def get_splitters(self) -> list[Splitter]:
        """Restituisce la lista degli splitter nel circuito."""
        return self.splitters

    def connect_blocks(self, source_block: 'Subcircuit', dest_block: 'Subcircuit',
                       source_output_names: list[str], dest_input_names: list[str]):
        """
        Connette i nodi di output di un blocco ai nodi di input di un altro blocco.
        Questo metodo assicura che i nodi siano gli stessi nel sistema MNA.
        Args:
            source_block (Subcircuit): Il blocco sorgente.
            dest_block (Subcircuit): Il blocco di destinazione.
            source_output_names (list[str]): Nomi dei nodi di output del blocco sorgente da connettere.
            dest_input_names (list[str]): Nomi dei nodi di input del blocco di destinazione da connettere.
        """
        if len(source_output_names) != len(dest_input_names):
            raise ValueError("Il numero di nodi di output sorgente deve corrispondere al numero di nodi di input destinazione.")

        for src_name, dest_name in zip(source_output_names, dest_input_names):
            # Assicurati che i nodi esistano in entrambi i sottocircuiti
            if src_name not in source_block.nodes:
                raise ValueError(f"Nodo '{src_name}' non trovato nel sottocircuito sorgente '{source_block.name}'.")
            if dest_name not in dest_block.nodes:
                raise ValueError(f"Nodo '{dest_name}' non trovato nel sottocircuito destinazione '{dest_block.name}'.")

            # Ottieni l'ID del nodo sorgente nel circuito principale
            src_node_id = self.nodes[source_block.nodes[src_name]]
            
            # Riassegna il nodo di destinazione all'ID del nodo sorgente
            # Questo unisce i due nodi nel circuito principale
            old_dest_node_id = self.nodes[dest_block.nodes[dest_name]]
            
            # Aggiorna il nodo di destinazione per tutti i componenti che lo usano
            for comp in self.components:
                # Se il componente ha nodi esplicitamente mappati per nome (come i pin_names)
                # o se i suoi node_ids sono una tupla che possiamo modificare
                if comp.node_ids is not None:
                    new_node_ids = list(comp.node_ids)
                    updated = False
                    for i, node_id in enumerate(new_node_ids):
                        if node_id == old_dest_node_id:
                            new_node_ids[i] = src_node_id
                            updated = True
                    if updated:
                        comp._set_node_ids(tuple(new_node_ids))
            
            # Aggiorna la mappatura principale del circuito
            # Rimuovi il vecchio nodo di destinazione e mappa il suo nome al nuovo ID
            del self.nodes[dest_block.nodes[dest_name]]
            self.nodes[dest_block.nodes[dest_name]] = src_node_id
            self.node_map[dest_block.nodes[dest_name]] = src_node_id

            print(f"Connesso '{source_block.name}:{src_name}' (ID {src_node_id}) a '{dest_block.name}:{dest_name}' (ex ID {old_dest_node_id}, ora {src_node_id}).")

