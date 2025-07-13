# circuit_modules/InputTransformerCircuit.py

from circuit_solver.circuit import Circuit
from components.inductor import Inductor
from components.resistor import Resistor

class InputTransformerCircuit(Circuit):
    """
    Circuito di ingresso con trasformatore.
    Modella un trasformatore ideale con avvolgimenti primari e secondari.
    Un trasformatore in MNA richiede tipicamente variabili ausiliarie per le correnti
    degli avvolgimenti e relazioni di tensione/corrente tra primario e secondario.
    Per la simulazione, semplificheremo con induttanze.
    """
    def __init__(self, name="InputTransformer", turns_ratio=1.0, primary_inductance_H=10.0, secondary_inductance_H=10.0, primary_resistance_ohm=10.0, secondary_resistance_ohm=10.0, sample_rate=48000):
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
        # Nodi del lato primario
        self.add_node("Input_P1") # Ingresso del primario
        self.add_node("Input_P2") # Altro capo del primario (spesso a GND o riferimento)

        # Nodi del lato secondario
        self.add_node("Output_S1") # Uscita del secondario
        self.add_node("Output_S2") # Altro capo del secondario (spesso a GND o riferimento)

        # GND è già gestito dalla classe base Circuit

    def _add_components(self):
        """Aggiunge i componenti specifici del trasformatore."""
        # Per un modello semplice di trasformatore:
        # Useremo due induttori (primario e secondario) e due resistori di avvolgimento.
        # L'accoppiamento magnetico verrà gestito implicitamente o necessiterà di equazioni MNA aggiuntive.
        # Questo è un modello semplificato, non un trasformatore accoppiato.

        # Resistenza dell'avvolgimento primario
        self.add_component(Resistor(self.primary_resistance), "Input_P1", "TP1_Internal")
        # Induttanza del primario
        self.add_component(Inductor(self.primary_inductance, sample_rate=self.sample_rate), "TP1_Internal", "Input_P2")

        # Resistenza dell'avvolgimento secondario
        self.add_component(Resistor(self.secondary_resistance), "Output_S1", "TS1_Internal")
        # Induttanza del secondario
        self.add_component(Inductor(self.secondary_inductance, sample_rate=self.sample_rate), "TS1_Internal", "Output_S2")

        # Nota: L'accoppiamento mutuo (M) tra gli induttori primario e secondario
        # è cruciale per un trasformatore reale. La MNA lo gestisce con equazioni
        # che coinvolgono le correnti e le derivate delle correnti degli induttori.
        # Queste equazioni dovranno essere aggiunte al solutore MNA o al componente Inductor
        # se si vuole un modello di trasformatore completo (generalmente, non si usa Inductor singolo).
        # Per ora, è una rappresentazione simbolica di "c'è un trasformatore".
        # Le sue equazioni MNA specifiche per il rapporto di spire e l'accoppiamento
        # dovrebbero essere gestite direttamente nel solver o in una classe Transformer dedicata.

    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        # Collego Input_P2 e Output_S2 a GND per un trasformatore single-ended
        # Se fosse bilanciato, Output_S2 andrebbe a un altro nodo.
        self.connect_nodes("Input_P2", "GND")
        self.connect_nodes("Output_S2", "GND")

    # Metodi per accedere ai nodi di input/output specifici del trasformatore
    def get_input_node(self):
        return self.get_node_id("Input_P1")

    def get_output_node(self):
        return self.get_node_id("Output_S1")
