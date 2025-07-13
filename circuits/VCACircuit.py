# circuit_modules/VCACircuit.py

from circuit_solver.circuit import Circuit
from components.mosfet import MOSFET # Useremo un MOSFET come elemento controllato

class VCACircuit(Circuit):
    """
    Circuito di un amplificatore controllato in tensione (VCA).
    Il guadagno è controllato da una tensione esterna.
    Modellato qui usando un MOSFET come resistore variabile in serie.
    """
    def __init__(self, name="VCA", mosfet_vt=1.0, mosfet_kn=0.001, series_resistance=1000.0):
        super().__init__(name)
        self.mosfet_vt = float(mosfet_vt)
        self.mosfet_kn = float(mosfet_kn)
        self.series_resistance = float(series_resistance)

        self.mosfet_component = None # Riferimento al componente MOSFET

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il VCA."""
        self.add_node("Audio_Input")
        self.add_node("Audio_Output")
        self.add_node("Control_Voltage_Input") # Il nodo di controllo (Gate del MOSFET)

        # Nodi interni per il MOSFET
        self.add_node("MOSFET_Drain")
        self.add_node("MOSFET_Source")
        self.add_node("MOSFET_Gate") # Questo è Control_Voltage_Input

    def _add_components(self):
        """Aggiunge i componenti del VCA."""
        # MOSFET come resistore variabile
        # Il Drain e Source sono in serie con il percorso audio.
        # Il Gate riceve la tensione di controllo.
        self.mosfet_component = MOSFET(Vt=self.mosfet_vt, Kn=self.mosfet_kn)
        self.add_component(self.mosfet_component, "MOSFET_Drain", "MOSFET_Source", "MOSFET_Gate")

        # Il MOSFET è in serie con il segnale audio (o in configurazione shunt).
        # Per un VCA, spesso si usa una configurazione a divisore di tensione.
        # Ad esempio: Ingresso -> R_series -> MOSFET -> Output -> R_shunt -> GND
        # dove il MOSFET è R_shunt o R_series.
        # Qui lo posizioniamo in serie. Avrai bisogno di un altro resistore per un divisore.
        
        # Resistor in serie al MOSFET per formare un divisore (se necessario)
        # self.add_component(Resistor(self.series_resistance), "Audio_Input", "MOSFET_Drain")
        # Connettiamo l'Output del VCA al Source del MOSFET.
        # Il controllo Gate è separato.

    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        # Collega il Gate del MOSFET al nodo di Control_Voltage_Input
        self.connect_nodes("MOSFET_Gate", "Control_Voltage_Input")
        
        # Connettiamo l'Input audio al Drain del MOSFET
        self.connect_nodes("Audio_Input", "MOSFET_Drain")
        # Connettiamo il Source del MOSFET all'Output audio
        self.connect_nodes("MOSFET_Source", "Audio_Output")
        # Emettiamo l'audio di uscita dal Source del MOSFET (configurazione source follower, o resistore)
        # o più tipicamente, il MOSFET è in un divisore di tensione.
        
        # Per un VCA come resistore variabile in serie:
        # Audio_Input --- R_fixed --- MOSFET_Drain -- MOSFET_Source --- Audio_Output --- (Load)
        # MOSFET_Gate si collega a Control_Voltage_Input
        # Qui il MOSFET è direttamente tra Audio_Input e Audio_Output, assumendo una configurazione
        # dove la sua resistenza Vds è l'elemento controllato.

        # Aggiungiamo un resistore di carico all'uscita per definire una tensione
        # self.add_component(Resistor(self.series_resistance), "Audio_Output", "GND") # Esempio di carico

    def set_control_voltage(self, control_voltage):
        """
        Imposta la tensione di controllo per il Gate del MOSFET.
        Questo metodo non aggiorna direttamente il MOSFET, ma è un'indicazione
        di dove la tensione di controllo verrebbe applicata nel sistema MNA.
        """
        # La tensione di controllo verrebbe applicata al nodo "Control_Voltage_Input"
        # nel solutore MNA, che a sua volta è collegato al Gate del MOSFET.
        # Non c'è bisogno di un'azione qui, dato che l'MNA lo gestirà.
        pass

    def get_audio_input_node(self):
        return self.get_node_id("Audio_Input")

    def get_audio_output_node(self):
        return self.get_node_id("Audio_Output")

    def get_control_voltage_node(self):
        return self.get_node_id("Control_Voltage_Input")

    def get_mosfet_component(self):
        return self.mosfet_component
