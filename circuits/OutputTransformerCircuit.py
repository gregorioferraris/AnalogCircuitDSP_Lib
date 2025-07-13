# circuit_modules/OutputTransformerCircuit.py

from circuit_solver.circuit import Circuit
from components.inductor import Inductor
from components.resistor import Resistor

class OutputTransformerCircuit(Circuit):
    """
    Circuito di uscita con trasformatore.
    Simile all'ingresso, ma la prospettiva è l'adattamento di impedenza verso un carico.
    """
    def __init__(self, name="OutputTransformer", turns_ratio=1.0, primary_inductance_H=10.0, secondary_inductance_H=10.0, primary_resistance_ohm=10.0, secondary_resistance_ohm=10.0, sample_rate=48000):
        super().__init__(name)
        self.turns_ratio = float(turns_ratio) # Rapporto spire Ns/Np
        self.primary_inductance = float(primary_inductance_H)
        self.secondary_inductance = float(secondary_inductance_H)
        self.primary_resistance = float(primary_resistance_ohm)
        self.secondary_resistance = float(secondary_resistance_ohm)
        self.sample_rate = sample_rate

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per questo circuito."""
        self.add_node("Input_P1") # Ingresso del primario (dalla fase di amplificazione)
        self.add_node("Input_P2") # Altro capo del primario (spesso a GND o alimentazione)

        self.add_node("Output_S1") # Uscita del secondario (verso il carico)
        self.add_node("Output_S2") # Altro capo del secondario (spesso a GND o riferimento)

    def _add_components(self):
        """Aggiunge i componenti specifici del trasformatore."""
        self.add_component(Resistor(self.primary_resistance), "Input_P1", "TP1_Internal")
        self.add_component(Inductor(self.primary_inductance, sample_rate=self.sample_rate), "TP1_Internal", "Input_P2")

        self.add_component(Resistor(self.secondary_resistance), "Output_S1", "TS1_Internal")
        self.add_component(Inductor(self.secondary_inductance, sample_rate=self.sample_rate), "TS1_Internal", "Output_S2")

        # Vedi note sul modello di trasformatore nel file InputTransformerCircuit.py

    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        self.connect_nodes("Input_P2", "GND") # Tipico per trasformatori in configurazione comune
        self.connect_nodes("Output_S2", "GND") # Output single-ended

    def get_input_node(self):
        return self.get_node_id("Input_P1")

    def get_output_node(self):
        return self.get_node_id("Output_S1")
