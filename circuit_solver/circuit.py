import numpy as np
from components.component import Component
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor
from components.voltage_source import VoltageSource
from components.current_source import CurrentSource
from components.diode import Diode
from components.transformer import Transformer
from components.splitter import Splitter
from circuit_solver.subcircuit import Subcircuit # Nuovo import

class Circuit:
    def __init__(self):
        self.components = []
        self.node_names = ['0'] # Nodo di terra per convenzione
        self.node_map = {'0': 0}
        self.num_nodes = 1
        self.voltage_source_count = 0 
        self.transformer_count = 0    
        self.splitter_current_count = 0 

    def add_component(self, component: Component):
        self.components.append(component)
        
        # Registra i nomi dei nodi e assegna ID numerici se non già presenti
        for node_name in component.node_names:
            if node_name not in self.node_map:
                self.node_map[node_name] = self.num_nodes
                self.node_names.append(node_name)
                self.num_nodes += 1
        
        # Assegna gli ID numerici ai nodi del componente
        component.node_ids = tuple(self.node_map[name] for name in component.node_names)

        # Pre-assegna gli indici per le correnti delle sorgenti di tensione
        if isinstance(component, VoltageSource):
            component._set_current_index(self.num_nodes + self.voltage_source_count)
            self.voltage_source_count += 1
        
        # Pre-assegna gli indici per le correnti dei trasformatori
        elif isinstance(component, Transformer):
            ip_idx = self.num_nodes + self.voltage_source_count + (self.transformer_count * 2)
            is_idx = ip_idx + 1
            component._set_current_indices(ip_idx, is_idx)
            self.transformer_count += 1
        
        # Pre-assegna gli indici per le correnti dello splitter
        elif isinstance(component, Splitter):
            start_idx = self.num_nodes + self.voltage_source_count + (self.transformer_count * 2) + self.splitter_current_count
            indices = [start_idx + i for i in range(component.num_outputs)]
            component._set_output_current_indices(indices)
            self.splitter_current_count += component.num_outputs


    def get_num_total_equations(self):
        """Restituisce il numero totale di equazioni (nodi + correnti Vs + correnti Trasformatore + correnti Splitter)."""
        return self.num_nodes + self.voltage_source_count + (self.transformer_count * 2) + self.splitter_current_count

    def connect_blocks(self, source_block: Subcircuit, dest_block: Subcircuit,
                       source_output_names: list[str] = None, dest_input_names: list[str] = None):
        """
        Connette i nodi di uscita di un blocco (Subcircuit) ai nodi di ingresso di un altro blocco.
        Questo metodo integra i componenti di entrambi i sottocircuiti nel circuito principale
        e unisce i nodi specificati.

        Args:
            source_block (Subcircuit): Il blocco sorgente (es. un preamplificatore).
            dest_block (Subcircuit): Il blocco di destinazione (es. un amplificatore di potenza).
            source_output_names (list[str], optional): Nomi specifici dei nodi di uscita del blocco sorgente da connettere.
                                                       Se None, usa tutti i nodi di uscita predefiniti del blocco sorgente.
            dest_input_names (list[str], optional): Nomi specifici dei nodi di ingresso del blocco di destinazione da connettere.
                                                    Se None, usa tutti i nodi di ingresso predefiniti del blocco di destinazione.
        Raises:
            ValueError: Se il numero di nodi da connettere non corrisponde o i nodi non sono definiti.
        """
        # Integra i componenti del blocco sorgente nel circuito principale
        for comp in source_block.components:
            self.add_component(comp)
        
        # Integra i componenti del blocco di destinazione nel circuito principale
        for comp in dest_block.components:
            self.add_component(comp)

        # Determina i nodi da connettere
        outputs_to_connect = source_output_names if source_output_names is not None else source_block.get_output_node_names()
        inputs_to_connect = dest_input_names if dest_input_names is not None else dest_block.get_input_node_names()

        if len(outputs_to_connect) != len(inputs_to_connect):
            raise ValueError(f"Il numero di nodi di uscita ({len(outputs_to_connect)}) del blocco sorgente '{source_block.name}' "
                             f"non corrisponde al numero di nodi di ingresso ({len(inputs_to_connect)}) del blocco di destinazione '{dest_block.name}'.")

        # Esegui le connessioni
        for i in range(len(outputs_to_connect)):
            src_node_name = outputs_to_connect[i]
            dest_node_name = inputs_to_connect[i]

            # Ottieni gli ID MNA dei nodi nel contesto del circuito principale
            # Se i nodi non sono stati aggiunti (dovrebbero esserlo tramite add_component), questo fallirà
            if src_node_name not in self.node_map:
                raise ValueError(f"Nodo di uscita '{src_node_name}' del blocco '{source_block.name}' non trovato nel circuito principale.")
            if dest_node_name not in self.node_map:
                raise ValueError(f"Nodo di ingresso '{dest_node_name}' del blocco '{dest_block.name}' non trovato nel circuito principale.")

            # Connetti i nodi: forza il nodo di destinazione ad avere lo stesso ID del nodo sorgente
            # Questo è il "cable" che unisce i pezzi
            # Utilizziamo wire_nodes per gestire la mappatura
            self.wire_nodes(src_node_name, dest_node_name)
            
            print(f"Connesso '{source_block.name}.{src_node_name}' a '{dest_block.name}.{dest_node_name}'.")

    def wire_nodes(self, *node_names_to_connect: str):
        """
        Connette un insieme di nodi tra loro, forzandoli ad essere lo stesso nodo logico.
        Elettricamente, significa che tutti i nodi passati avranno lo stesso ID numerico.
        Questo metodo è robusto e gestisce il remapping dei componenti esistenti.
        """
        if not node_names_to_connect:
            return

        # Il primo nodo nel gruppo sarà il "master"
        master_node_name = node_names_to_connect[0]
        if master_node_name not in self.node_map:
            self.node_map[master_node_name] = self.num_nodes
            self.node_names.append(master_node_name)
            self.num_nodes += 1
        master_node_id = self.node_map[master_node_name]

        # Raccogli i nodi che devono essere remappati
        nodes_to_remap = []
        for i in range(1, len(node_names_to_connect)):
            current_node_name = node_names_to_connect[i]
            if current_node_name in self.node_map and self.node_map[current_node_name] != master_node_id:
                nodes_to_remap.append((current_node_name, self.node_map[current_node_name]))
            self.node_map[current_node_name] = master_node_id
            if current_node_name not in self.node_names:
                self.node_names.append(current_node_name)

        # Esegui il remapping degli ID nei componenti già aggiunti
        if nodes_to_remap:
            print(f"Eseguendo remapping per nodi: {nodes_to_remap} al master ID {master_node_id}")
            for comp in self.components:
                # node_ids è una tupla, quindi la convertiamo in lista per modificarla
                new_node_ids = list(comp.node_ids) 
                
                changed = False
                for i, node_id in enumerate(new_node_ids):
                    for old_name, old_id in nodes_to_remap:
                        if node_id == old_id:
                            new_node_ids[i] = master_node_id
                            changed = True
                            break # Esci dal loop dei nodi da remappare per questo pin
                
                if changed:
                    comp.node_ids = tuple(new_node_ids)
                    # print(f"  Aggiornato {comp.name}: nuovi node_ids {comp.node_ids}")

        print(f"Nodi {node_names_to_connect} connessi al nodo logico '{master_node_name}' (ID: {master_node_id}).")

