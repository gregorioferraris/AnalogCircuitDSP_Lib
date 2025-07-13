# circuits/TubePowerAmpCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.pentode import Pentode # Il TUO componente Pentode aggiornato
# Potrebbe essere necessario un componente trasformatore ideale per l'accoppiamento di uscita
# from components.transformer import IdealTransformer # Se ne hai uno

import numpy as np

class TubePowerAmpCircuit(Circuit):
    """
    Circuito di un amplificatore di potenza valvolare Single-Ended semplificato.
    Utilizza un Pentodo con un carico resistivo (idealizzato) che simula l'impedenza del trasformatore
    o un trasformatore di uscita ideale per l'accoppiamento al diffusore.
    """
    def __init__(self, name="Tube_Power_Amp",
                 pentode_params=None, # Parametri del Pentodo
                 R_grid_input=100000.0, # Resistenza di ingresso griglia (Grid Leak Resistor)
                 R_cathode=500.0,    # Resistenza di catodo (Cathode Resistor for bias)
                 C_cathode=100.0e-6, # Condensatore di bypass catodo
                 R_screen_grid=1000.0, # Resistenza della griglia schermo (Screen Grid Resistor)
                 C_screen_grid=10.0e-6, # Condensatore di bypass della griglia schermo
                 R_anode_load=5000.0, # Carico ideale della placca (simulazione di impedenza del trasformatore)
                 C_coupling_in=0.1e-6, # Condensatore di accoppiamento in ingresso
                 V_power_supply_anode=350.0, # Tensione di alimentazione B+ per l'anodo
                 V_power_supply_screen=300.0, # Tensione di alimentazione per la griglia schermo (G2)
                 sample_rate=48000):
        super().__init__(name)
        self.sample_rate = sample_rate

        self.pentode_params = pentode_params if pentode_params is not None else {}
        self.R_grid_input_val = float(R_grid_input)
        self.R_cathode_val = float(R_cathode)
        self.C_cathode_val = float(C_cathode)
        self.R_screen_grid_val = float(R_screen_grid)
        self.C_screen_grid_val = float(C_screen_grid)
        self.R_anode_load_val = float(R_anode_load)
        self.C_coupling_in_val = float(C_coupling_in)
        self.V_power_supply_anode = float(V_power_supply_anode)
        self.V_power_supply_screen = float(V_power_supply_screen)

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per l'amplificatore di potenza."""
        self.add_node("Input")
        self.add_node("Output") # Nodo per l'uscita audio
        self.add_node("V_Power_Supply_Anode") # Nodo per l'alimentazione B+ dell'anodo
        self.add_node("V_Power_Supply_Screen") # Nodo per l'alimentazione della griglia schermo

        # Nodi della Valvola Pentodo
        self.add_node("Pentode_Anode")
        self.add_node("Pentode_Grid1") # Griglia di controllo
        self.add_node("Pentode_Grid2") # Griglia schermo
        self.add_node("Pentode_Grid3") # Griglia soppressore (normalmente a catodo/GND)
        self.add_node("Pentode_Cathode")

        # Nodi intermedi
        self.add_node("Input_Coupling_Node")
        self.add_node("Screen_Grid_Point") # Nodo tra R_screen_grid e C_screen_grid

    def _add_components(self):
        """Aggiunge i componenti (Pentodo, R, C, Carico)."""
        # Sorgenti di alimentazione B+
        self.add_voltage_source(self.V_power_supply_anode, "V_Power_Supply_Anode", "GND", name="B_Plus_Supply_Anode")
        self.add_voltage_source(self.V_power_supply_screen, "V_Power_Supply_Screen", "GND", name="Screen_Supply")

        # Valvola Pentodo
        self.pentode = Pentode(name="PowerAmp_Pentode", **self.pentode_params)
        self.pentode.set_nodes(anode_node_id=self.get_node_id("Pentode_Anode"),
                               grid1_node_id=self.get_node_id("Pentode_Grid1"),
                               grid2_node_id=self.get_node_id("Pentode_Grid2"),
                               cathode_node_id=self.get_node_id("Pentode_Cathode"))
        self.nonlinear_components.append(self.pentode)

        # Componenti di ingresso (Griglia di controllo G1)
        self.add_component(Capacitor(self.C_coupling_in_val, sample_rate=self.sample_rate),
                           "Input", "Input_Coupling_Node", name="C_coupling_in_PA")
        self.add_component(Resistor(self.R_grid_input_val),
                           "Input_Coupling_Node", "GND", name="R_grid_input_PA")
        
        # Connessione della griglia di controllo (G1)
        self.connect_nodes("Pentode_Grid1", "Input_Coupling_Node")

        # Componenti del catodo (bias e bypass)
        self.add_component(Resistor(self.R_cathode_val),
                           "Pentode_Cathode", "GND", name="R_cathode_PA")
        self.add_component(Capacitor(self.C_cathode_val, sample_rate=self.sample_rate),
                           "Pentode_Cathode", "GND", name="C_cathode_PA")

        # Componenti della griglia schermo (G2)
        self.add_component(Resistor(self.R_screen_grid_val),
                           "V_Power_Supply_Screen", "Screen_Grid_Point", name="R_screen_grid")
        self.add_component(Capacitor(self.C_screen_grid_val, sample_rate=self.sample_rate),
                           "Screen_Grid_Point", "GND", name="C_screen_grid")
        
        # Connessione della griglia schermo (G2)
        self.connect_nodes("Pentode_Grid2", "Screen_Grid_Point")

        # Connessione della griglia soppressore (G3) a catodo (o GND)
        self.connect_nodes("Pentode_Grid3", "Pentode_Cathode") # o "GND" se il catodo è a GND

        # Carico di placca (simula l'impedenza del trasformatore di uscita)
        self.add_component(Resistor(self.R_anode_load_val),
                           "V_Power_Supply_Anode", "Pentode_Anode", name="R_anode_load_PA")

        # L'uscita è prelevata dalla placca. In un sistema reale, andrebbe a un trasformatore.
        # Per semplicità, la placca è il nodo di uscita, assumendo un blocco DC esterno.

    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        self.connect_nodes("Output", "Pentode_Anode") # Uscita direttamente dalla placca

    def get_input_node(self):
        return self.get_node_id("Input")

    def get_output_node(self):
        return self.get_node_id("Output")
