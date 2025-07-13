# circuit_modules/VariMuCircuit.py

from circuit_solver.circuit import Circuit
from components.pentode import Pentode # I pentodi sono spesso usati per Vari-Mu
from components.triode import Triode # A volte anche triodi speciali o configurazioni
from components.resistor import Resistor
from components.capacitor import Capacitor

class VariMuCircuit(Circuit):
    """
    Circuito di un compressore Vari-Mu (Variable-Mu).
    Utilizza la caratteristica di transconduttanza variabile delle valvole (spesso pentodi)
    per ottenere la compressione del guadagno. La tensione di controllo (bias della griglia)
    modifica il guadagno della valvola.
    """
    def __init__(self, name="VariMuCompressorStage", tube_type="pentode", operating_point_vbias=-5.0, anode_resistor=22000.0, cathode_resistor=1000.0, coupling_capacitor=1.0e-6, sample_rate=48000):
        super().__init__(name)
        self.tube_type = tube_type.lower()
        self.operating_point_vbias = float(operating_point_vbias) # Tensione di bias della griglia
        self.anode_resistor = float(anode_resistor)
        self.cathode_resistor = float(cathode_resistor)
        self.coupling_capacitor = float(coupling_capacitor)
        self.sample_rate = sample_rate

        self.gain_tube = None # Riferimento alla valvola di guadagno
        self.cathode_capacitor_value = 100e-6 # Valore tipico per bypass catodico

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per lo stadio Vari-Mu."""
        self.add_node("Audio_Input")
        self.add_node("Audio_Output")
        self.add_node("Control_Voltage_Input") # Per la tensione di bias variabile sulla griglia

        # Nodi interni della valvola
        self.add_node("Tube_Anode")
        self.add_node("Tube_Grid1") # Griglia di controllo
        self.add_node("Tube_Grid2") # Griglia schermo (per pentodi)
        self.add_node("Tube_Cathode")

        # Nodi per l'alimentazione (B+)
        self.add_node("B_PLUS") # Tensione positiva di alimentazione (es. 250V)


    def _add_components(self):
        """Aggiunge i componenti dello stadio Vari-Mu."""
        # Resistor di accoppiamento in ingresso per il bias
        self.add_component(Resistor(self.anode_resistor), "B_PLUS", "Tube_Anode")

        # La valvola principale di guadagno (Pentodo o Triodo)
        if self.tube_type == "pentode":
            # Parametri di esempio per un pentodo Vari-Mu (es. 6BA6, 6386)
            # Questi valori sono molto generici e necessiterebbero calibrazione reale
            self.gain_tube = Pentode(mu=20.0, Kp=0.2, X=1.5, Kg1=2.0, Kg2=10.0)
            self.add_component(self.gain_tube, "Tube_Anode", "Tube_Grid1", "Tube_Grid2", "Tube_Cathode")
            # Resistor e Capacitor per Griglia 2 (Screen Grid) - tipicamente a un potenziale fisso
            self.add_component(Resistor(100000.0), "B_PLUS", "Tube_Grid2") # Resistor di limitazione per G2
            self.add_component(Capacitor(1.0e-7, sample_rate=self.sample_rate), "Tube_Grid2", "GND") # Bypass di G2

        elif self.tube_type == "triode":
            # Parametri di esempio per un triodo Vari-Mu (meno comune per alta compressione)
            self.gain_tube = Triode(mu=20.0, Kp=0.5, X=1.5, Kg1=5.0)
            self.add_component(self.gain_tube, "Tube_Anode", "Tube_Grid1", "Tube_Cathode")
            # Per un triodo, non c'è Tube_Grid2, quindi le connessioni saranno diverse.

        else:
            raise ValueError(f"Tipo di valvola '{self.tube_type}' non supportato per Vari-Mu. Usa 'pentode' o 'triode'.")

        # Resistor di catodo per auto-bias (o bias fisso se non bypassato)
        self.add_component(Resistor(self.cathode_resistor), "Tube_Cathode", "GND")
        # Condensatore di bypass del catodo (per guadagno AC)
        self.add_component(Capacitor(self.cathode_capacitor_value, sample_rate=self.sample_rate), "Tube_Cathode", "GND")

        # Condensatore di accoppiamento in ingresso (per bloccare il DC bias)
        self.add_component(Capacitor(self.coupling_capacitor, sample_rate=self.sample_rate), "Audio_Input", "Tube_Grid1")

        # Condensatore di accoppiamento in uscita (per bloccare il DC bias)
        self.add_component(Capacitor(self.coupling_capacitor, sample_rate=self.sample_rate), "Tube_Anode", "Audio_Output")


    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        # Il segnale di controllo (bias DC) viene sommato al bias fisso o direttamente applicato alla griglia
        # Qui, "Control_Voltage_Input" è il nodo a cui si collega il segnale di controllo (es. da un raddrizzatore).
        # Assumiamo che ci sia un resistore di griglia (Grid Leak Resistor) per impostare il bias.
        self.add_component(Resistor(470000.0), "Tube_Grid1", "Control_Voltage_Input") # Esempio di resistore di griglia


    def set_control_voltage(self, control_voltage):
        """
        Imposta la tensione di controllo sulla griglia della valvola.
        Questo metodo non aggiorna direttamente la valvola, ma è un'indicazione
        di dove la tensione di controllo verrebbe applicata nel sistema MNA.
        """
        # La tensione di controllo (che varierà nel tempo in base al segnale audio)
        # verrà applicata al nodo "Control_Voltage_Input" nel solutore MNA.
        # Questa tensione modifica il bias sulla griglia e quindi la transconduttanza della valvola.
        pass

    def get_audio_input_node(self):
        return self.get_node_id("Audio_Input")

    def get_audio_output_node(self):
        return self.get_node_id("Audio_Output")

    def get_control_voltage_node(self):
        return self.get_node_id("Control_Voltage_Input")

    def get_gain_tube(self):
        return self.gain_tube
