# circuits/LowPassRCFilterCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor

class LowPassRCFilterCircuit(Circuit):
    """
    Circuito di un filtro passa-basso RC del primo ordine.
    Composto da un resistore e un condensatore.
    """
    def __init__(self, name="RC_LowPass", resistance_ohm=1000.0, capacitance_farad=1.0e-7, sample_rate=48000):
        super().__init__(name)
        self.R = float(resistance_ohm)
        self.C = float(capacitance_farad)
        self.sample_rate = sample_rate

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il filtro passa-basso RC."""
        self.add_node("Input")
        self.add_node("Output")
        # GND è già gestito

    def _add_components(self):
        """Aggiunge i componenti R e C."""
        self.add_component(Resistor(self.R), "Input", "Internal_Node")
        self.add_component(Capacitor(self.C, sample_rate=self.sample_rate), "Internal_Node", "GND")
        
        # L'output è preso dal condensatore
        # Non aggiungiamo qui un componente, ma usiamo l'Internal_Node come Output
        self.connect_nodes("Output", "Internal_Node")

    def _connect_nodes(self):
        """Non sono necessarie connessioni aggiuntive oltre quelle definite in add_components."""
        pass

    def get_input_node(self):
        return self.get_node_id("Input")

    def get_output_node(self):
        return self.get_node_id("Output")

    def get_cutoff_frequency(self):
        """Calcola la frequenza di taglio (-3dB) del filtro."""
        return 1.0 / (2.0 * np.pi * self.R * self.C)
