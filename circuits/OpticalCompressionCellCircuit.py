# circuit_modules/OpticalCompressionCellCircuit.py

from circuit_solver.circuit import Circuit
from components.ldr import LDR
from components.led import LED
from components.resistor import Resistor

class OpticalCompressionCellCircuit(Circuit):
    """
    Circuito della cellula ottica di compressione (Vactrol-like).
    Comprende una Light Dependent Resistor (LDR) e un LED accoppiati otticamente.
    Il LED è pilotato da una tensione di controllo (audio rettificato o DC),
    e la LDR varia la sua resistenza in base alla luce emessa dal LED, influenzando
    il percorso del segnale audio.
    """
    def __init__(self, name="OpticalCompressionCell", ldr_r_dark=1e6, ldr_r_light=100, led_saturation_current=1e-20, led_emission_coefficient=3.4, led_current_limiting_resistor=1000.0, sample_rate=48000):
        super().__init__(name)
        self.ldr_r_dark = float(ldr_r_dark)
        self.ldr_r_light = float(ldr_r_light)
        self.led_is = float(led_saturation_current)
        self.led_n = float(led_emission_coefficient)
        self.led_resistor = float(led_current_limiting_resistor)
        self.sample_rate = sample_rate

        self.ldr_component = None # Riferimento al componente LDR per aggiornare il livello di luce
        self.led_component = None # Riferimento al componente LED

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per la cellula ottica."""
        # Nodi per il LED
        self.add_node("LED_Control_Input") # Tensione di controllo per il LED
        self.add_node("LED_Anode")
        self.add_node("LED_Cathode")

        # Nodi per la LDR (parte del percorso audio)
        self.add_node("LDR_Input") # Ingresso audio alla LDR
        self.add_node("LDR_Output") # Uscita audio dalla LDR

        # GND è già gestito

    def _add_components(self):
        """Aggiunge i componenti della cellula ottica."""
        # Resistor in serie al LED (per limitare la corrente e renderlo controllabile)
        self.add_component(Resistor(self.led_resistor), "LED_Control_Input", "LED_Anode")
        # LED stesso
        self.led_component = LED(saturation_current=self.led_is, emission_coefficient=self.led_n)
        self.add_component(self.led_component, "LED_Anode", "LED_Cathode")

        # LDR
        self.ldr_component = LDR(R_dark=self.ldr_r_dark, R_light=self.ldr_r_light, sample_rate=self.sample_rate)
        self.add_component(self.ldr_component, "LDR_Input", "LDR_Output")


    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        # Il catodo del LED è solitamente a GND
        self.connect_nodes("LED_Cathode", "GND")

        # LDR opera in serie o in shunt nel percorso audio.
        # Non la colleghiamo internamente qui a un percorso audio completo,
        # ma forniamo i nodi LDR_Input e LDR_Output per l'integrazione nel compressore.

    def set_control_voltage(self, control_voltage):
        """
        Imposta la tensione di controllo per il LED e, di conseguenza, il livello di luce.
        Questo metodo verrà chiamato dal circuito di controllo del compressore.
        """
        # Questo è un placeholder. Idealmente, il solutore MNA risolverà la corrente del LED
        # data la tensione su LED_Control_Input.
        # Per un modello disaccoppiato, potremmo calcolare la corrente del LED e poi la luce.

        # Qui, usiamo una semplificazione: simuliamo la corrente del LED data la tensione
        # ai suoi capi (che il solutore calcolerà).
        # Per ora, si può calcolare la corrente del LED data la tensione e la resistenza in serie.
        # Questa corrente determina il livello di luce per la LDR.

        # Una stima semplificata della corrente del LED:
        # Se V_control > V_forward_LED, allora I_LED = (V_control - V_forward_LED) / R_series
        # Altrimenti I_LED = 0
        
        # Questa logica di "controllo" deve essere gestita all'esterno dal circuito principale.
        # Il solutore MNA calcola le correnti e le tensioni.
        # La LDR deve essere aggiornata con il "livello di luce" che dipende dalla corrente del LED.
        
        # Una semplificazione: il "livello di luce" per la LDR è una funzione della tensione di controllo.
        # light_level = max(0, (control_voltage - self.led_component.Vf_typical) / self.led_resistor) # Corrente stimata del LED
        
        # O più semplicemente, la luce è proporzionale alla tensione di controllo (sotto una certa soglia)
        # Una mappatura non lineare è più realistica
        max_light_level = 1000.0 # Valore arbitrario massimo per la luce
        # La luce emessa è proporzionale alla corrente del LED
        # I_led_current = self.led_component.calculate_current(control_voltage - self.led_component.Vf_typical) # NON CORRETTO, LED_Control_Input è un nodo.
        
        # La corrente che attraversa il LED si risolve nel sistema MNA.
        # Ma per aggiornare la LDR, abbiamo bisogno del "livello di luce".
        # Qui potremmo creare una mappatura diretta dalla tensione di controllo al livello di luce,
        # oppure calcolare la corrente del LED dal sistema MNA risolto e usarla per l'LDR.
        
        # Per il momento, assumiamo che 'light_level' sia una funzione proporzionale alla tensione di controllo,
        # magari saturata, e passiamo questa alla LDR.
        # Questo è un HACK, la corrente del LED dovrebbe venire dalla soluzione MNA.
        # Supponiamo un mapping semplice da control_voltage a light_level.
        # Un valore di esempio: 0V -> 0 luce, 5V -> 1000 luce
        
        # Mappatura approssimativa:
        light_level = max(0, control_voltage * (max_light_level / 5.0)) # Scala 0-5V a 0-1000
        
        # Aggiorna la resistenza della LDR in base a questo livello di luce "esterno"
        self.ldr_component.get_resistance(light_level) # get_resistance aggiorna la current_resistance interna


    # Metodi per accedere ai nodi di input/output specifici
    def get_ldr_input_node(self):
        return self.get_node_id("LDR_Input")

    def get_ldr_output_node(self):
        return self.get_node_id("LDR_Output")

    def get_led_control_node(self):
        return self.get_node_id("LED_Control_Input")

    def get_ldr_component(self):
        return self.ldr_component

    def get_led_component(self):
        return self.led_component
