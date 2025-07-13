# circuit_solver/circuit.py

import numpy as np

# Importa tutti i tipi di componenti che possono essere aggiunti al circuito
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor
from components.diode import Diode
from components.led import LED
from components.ldr import LDR # LDR è un po' speciale, è una resistenza variabile
from components.jfet import JFET
from components.bjt import BJT
from components.schottky_diode import SchottkyDiode
from components.zener_diode import ZenerDiode
from components.mosfet import MOSFET
from components.triode import Triode
from components.pentode import Pentode
from components.rectifier_tube import RectifierTube

# Potresti anche voler importare altre utilità se necessarie
# from utils.constants import DEFAULT_SAMPLE_RATE

class Circuit:
    """
    Rappresenta la struttura di un circuito, inclusi nodi e componenti.
    Questa classe si occupa di "assemblare" il circuito.
    """
    def __init__(self, name="Unnamed Circuit"):
        self.name = name
        self.nodes = {}             # Dizionario per memorizzare i nomi dei nodi e i loro indici numerici
        self.node_count = 0         # Contatore per assegnare indici unici ai nodi
        self.components = []        # Lista degli oggetti componente
        self.component_id_counter = 0 # Contatore per ID unici dei componenti

        # Nodi speciali
        self.ground_node_id = 0     # Il nodo 0 è sempre la massa (ground)
        self.add_node("GND")        # Aggiungi il nodo di massa inizialmente

        print(f"Circuito '{self.name}' inizializzato. GND è il nodo {self.ground_node_id}.")

    def add_node(self, node_name):
        """
        Aggiunge un nodo al circuito se non esiste già.
        Restituisce l'ID numerico del nodo.
        """
        if node_name not in self.nodes:
            self.nodes[node_name] = self.node_count
            self.node_count += 1
            # print(f"Nodo '{node_name}' aggiunto con ID: {self.nodes[node_name]}")
        return self.nodes[node_name]

    def get_node_id(self, node_name):
        """
        Restituisce l'ID numerico di un nodo dato il suo nome.
        Se il nodo non esiste, lo crea.
        """
        return self.add_node(node_name)

    def get_node_name(self, node_id):
        """
        Restituisce il nome di un nodo dato il suo ID numerico.
        """
        for name, id_val in self.nodes.items():
            if id_val == node_id:
                return name
        raise ValueError(f"Nodo con ID {node_id} non trovato.")

    def add_component(self, component_object, *node_names):
        """
        Aggiunge un componente al circuito e lo connette ai nodi specificati.
        Args:
            component_object: Un'istanza di una classe componente (es. Resistor(100)).
            *node_names: Nomi dei nodi a cui il componente è connesso.
                         Il numero di nodi dipende dal tipo di componente.
                         Es: Resistor('nodeA', 'nodeB'), Triode('grid', 'anode', 'cathode').
        """
        # Assegna un ID univoco al componente
        component_id = self.component_id_counter
        self.component_id_counter += 1
        setattr(component_object, 'component_id', component_id) # Aggiunge l'ID all'oggetto

        # Mappa i nomi dei nodi agli ID numerici
        connected_node_ids = tuple(self.get_node_id(name) for name in node_names)
        setattr(component_object, 'connected_nodes', connected_node_ids) # Associa i nodi

        self.components.append(component_object)
        print(f"Componente '{component_object.__class__.__name__}' (ID: {component_id}) aggiunto, connesso a nodi: {node_names} (IDs: {connected_node_ids}).")

        # Logica per gestire componenti speciali come LDR che sono resistenze variabili
        # Potresti voler salvare un riferimento diretto alla LDR per controllarla esternamente
        if isinstance(component_object, LDR):
            self.ldr_component = component_object
            print("LDR rilevata e memorizzata per il controllo esterno.")

    def get_num_nodes(self):
        """Restituisce il numero totale di nodi nel circuito (incluso GND)."""
        return self.node_count

    def get_components(self):
        """Restituisce la lista di tutti i componenti nel circuito."""
        return self.components

    def get_ground_node_id(self):
        """Restituisce l'ID del nodo di massa."""
        return self.ground_node_id

    def __str__(self):
        s = f"Circuito: {self.name}\n"
        s += f"Nodi ({self.node_count}): {self.nodes}\n"
        s += "Componenti:\n"
        for comp in self.components:
            s += f"  - {comp} (ID: {comp.component_id}, Connesso a: {comp.connected_nodes})\n"
        return s

# Esempio di utilizzo (per testare la struttura)
if __name__ == "__main__":
    my_circuit = Circuit("Semplice Amplificatore RC")

    # Aggiungi nodi (anche se add_component li aggiunge automaticamente se non esistono)
    my_circuit.add_node("Input")
    my_circuit.add_node("Output")
    my_circuit.add_node("VCC")

    # Aggiungi componenti
    # Resistore tra Input e VCC
    my_circuit.add_component(Resistor(10e3), "Input", "VCC")
    # Condensatore tra Input e GND
    my_circuit.add_component(Capacitor(1e-6), "Input", "GND")
    # Diodo tra Input e Output
    my_circuit.add_component(Diode(), "Input", "Output")
    # JFET
    my_circuit.add_component(JFET(), "DrainNode", "GateNode", "SourceNode")
    # Triodo (Anode, Grid, Cathode)
    my_circuit.add_component(Triode(), "AnodeNode", "GridNode", "CathodeNode")


    # Esempio di LDR che dovrà essere controllata esternamente
    my_circuit.add_component(LDR(), "LDR_Node1", "LDR_Node2")


    print("\n--- Stato del Circuito ---")
    print(my_circuit)

    # Verifica se la LDR è stata memorizzata
    if hasattr(my_circuit, 'ldr_component'):
        print(f"\nRiferimento LDR trovato: {my_circuit.ldr_component}")

    # Esempio di accesso ai nodi e componenti
    print(f"\nID del nodo 'Output': {my_circuit.get_node_id('Output')}")
    print(f"Nome del nodo con ID 1: {my_circuit.get_node_name(1)}")

    # Un'istruzione che fallirebbe per un nodo non esistente
    try:
        my_circuit.get_node_name(99)
    except ValueError as e:
        print(f"Errore previsto: {e}")
