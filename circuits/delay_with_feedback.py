# circuits/delay_with_feedback.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from circuits.analog_buffer import AnalogBuffer # Il buffer che abbiamo appena creato
from components.attenuator import Attenuator # L'attenuatore che abbiamo appena creato
from circuit_solver.circuit import Circuit
# Per la classe DelayLine, assumiamo che tu l'abbia come componente funzionale nel solver
# o che il solver abbia un modo per accedere al suo output e gestirlo.

class DelayWithFeedback(Circuit):
    def __init__(self, name="Delay_With_Feedback_Loop"):
        super().__init__(name)

        # --- Parametri del circuito ---
        # Parametri sorgente AC d'ingresso
        self.input_amplitude = 1.0 # V
        self.input_frequency = 1000 # Hz

        # Parametri del delay funzionale (gestito esternamente, ma i nodi sono qui)
        # Nota: Questi sono solo placeholder per i nodi, la logica del delay_line è nel solver
        # o in una classe helper collegata al solver.
        self.delay_time_seconds = 0.5 # 500 ms di ritardo iniziale
        self.feedback_attenuation_factor = 0.7 # 70% di feedback (da 0 a 1)

        # Condensatori di accoppiamento
        C_in_val = 1.0e-6
        C_out_val = 1.0e-6

        # Resistenza per la somma del feedback (mixer passivo)
        R_sum_val = 10e3

        # --- Nodi del circuito ---
        self.add_node('main_in')       # Ingresso principale del delay
        self.add_node('delay_input')   # Ingresso effettivo del modulo di ritardo (dopo la somma)
        self.add_node('delay_output')  # Uscita del modulo di ritardo
        self.add_node('feedback_loop_node') # Nodo all'interno del loop di feedback
        self.add_node('main_out')      # Uscita finale del circuito delay
        self.add_node('V_in_AC')       # Nodo per la sorgente di segnale AC

        # --- Componenti ---

        # 1. Sorgente di segnale d'ingresso AC
        self.add_component(VoltageSource(name="V_signal_input", nodes={'pos': 'V_in_AC', 'neg': 'gnd'}, voltage=0.0))
        self.add_component(Resistor(name="R_in_iso", node1='V_in_AC', node2='main_in', resistance=1e3)) # Resistenza isolante

        # 2. Buffer Analogico (prima dell'ingresso del delay per isolamento)
        # Qui useremo il tuo AnalogBuffer come sottocircuito.
        # Devi aggiungere i suoi componenti al circuito principale e connettere i nodi.
        self.buffer = AnalogBuffer(name="Input_Buffer")
        # Connetti l'ingresso del buffer al main_in del delay
        self.buffer_input_node = self.buffer.get_component("C_in").node1 # L'input effettivo del buffer è V_in_AC (prima di C_in)
        self.add_component(Capacitor(name="C_main_in_buffer", node1='main_in', node2=self.buffer.get_component("C_in").node1, capacitance=C_in_val))
        self.add_component(Resistor(name="R_main_in_buffer", node1='main_in', node2=self.buffer.get_component("RG").node1, resistance=1e3))
        
        # Aggiungi i componenti interni del buffer al circuito principale
        for comp in self.buffer.components:
            self.add_component(comp)
        
        # Rinominare il nodo di uscita del buffer per chiarezza nel contesto di questo circuito
        self.map_node(self.buffer.get_component("C_out").node1, 'buffered_input') # L'uscita del buffer è ora 'buffered_input'
        
        # 3. Stadio di Somma (Mixer Passivo per Input + Feedback)
        # Il segnale di ingresso originale (da 'buffered_input') si somma con il feedback ('feedback_loop_node').
        self.add_component(Resistor(name="R_sum_input", node1='buffered_input', node2='delay_input', resistance=R_sum_val))
        self.add_component(Resistor(name="R_sum_feedback", node1='feedback_loop_node', node2='delay_input', resistance=R_sum_val))
        self.add_component(Resistor(name="R_sum_to_gnd", node1='delay_input', node2='gnd', resistance=R_sum_val)) # Per stabilità


        # 4. Modulo di Ritardo (DelayLine - *simulato come un blocco funzionale dal solver*)
        # Qui NON aggiungiamo una componente fisica, ma definiamo i nodi di input/output che il solver userà.
        # Il solver prenderà V('delay_input') e fornirà V('delay_output') dopo un tempo delay_time_seconds.
        # Per la simulazione DC, questo blocco non fa nulla, è per la transitoria.
        # In un'implementazione reale, il solver dovrebbe avere accesso a questa 'DelayLine'
        # e usarla nel suo loop di integrazione.
        # Ad esempio: V_delay_out = self.delay_line_instance.update(V_delay_input, dt)
        # E poi impostare una sorgente di tensione V_delay_output_source = V_delay_out.
        
        # Per la simulazione in MNA, possiamo rappresentarlo con una sorgente di tensione
        # il cui valore è controllato esternamente dal nostro simulatore.
        # Questa sarà la "uscita" del delay che il solver deve calcolare dinamicamente.
        self.add_component(VoltageSource(name="V_delay_output_source", nodes={'pos': 'delay_output', 'neg': 'gnd'}, voltage=0.0))


        # 5. Buffer Analogico (all'uscita del delay, prima del feedback e dell'uscita principale)
        self.output_buffer = AnalogBuffer(name="Output_Buffer")
        self.add_component(Capacitor(name="C_delay_out_buffer", node1='delay_output', node2=self.output_buffer.get_component("C_in").node1, capacitance=C_in_val))
        self.add_component(Resistor(name="R_delay_out_buffer", node1='delay_output', node2=self.output_buffer.get_component("RG").node1, resistance=1e3))
        
        for comp in self.output_buffer.components:
            self.add_component(comp)
        
        # Rinominare il nodo di uscita del buffer per chiarezza
        self.map_node(self.output_buffer.get_component("C_out").node1, 'buffered_delay_out') # L'uscita del buffer è 'buffered_delay_out'


        # 6. Attenuatore (per il controllo del Feedback)
        self.feedback_attenuator = Attenuator(name="Feedback_Attenuator", 
                                              node1='buffered_delay_out', # Ingresso dall'uscita del delay (bufferizzata)
                                              node2='gnd', 
                                              node_tap='feedback_loop_node', # Uscita che va al nodo di somma
                                              attenuation_factor=self.feedback_attenuation_factor)
        for comp in self.feedback_attenuator.components:
            self.add_component(comp)

        # 7. Uscita Finale (potrebbe essere un altro mixer per Dry/Wet, ma qui solo il Wet delayed)
        self.add_component(Resistor(name="R_final_output", node1='buffered_delay_out', node2='main_out', resistance=1e3))
        self.add_component(Resistor(name="R_main_out_load", node1='main_out', node2='gnd', resistance=10e3))


# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    delay_circuit = DelayWithFeedback()
    
    print(f"Circuito: {delay_circuit.name}")
    print("Nodi:", delay_circuit.get_node_names())
    print("Componenti:")
    for comp in delay_circuit.components:
        print(f"  - {comp.name}: {type(comp).__name__} ({comp.name})")

    solver = MnaSolver(delay_circuit)

    print("\nSimulazione DC del Delay With Feedback Circuit...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in delay_circuit.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Verifica che i nodi AC siano a 0V in DC
        expected_zero_nodes = ['main_in', 'delay_input', 'delay_output', 'feedback_loop_node', 'main_out']
        for node_name in expected_zero_nodes:
            v_node_dc = node_voltages[solver._get_node_index(node_name)]
            if abs(v_node_dc) > 1e-6:
                print(f"ATTENZIONE: V({node_name}) DC non è a zero: {v_node_dc:.3f} V")
            else:
                print(f"V({node_name}) DC: {v_node_dc:.3f} V (ok)")

        # Verifica il funzionamento dell'attenuatore in DC (se V_signal_input è a 0, anche l'input dell'attenuatore è 0)
        v_buffered_delay_out_dc = node_voltages[solver._get_node_index('buffered_delay_out')]
        v_feedback_loop_node_dc = node_voltages[solver._get_node_index('feedback_loop_node')]
        print(f"\nVerifica Attenuatore (in DC, tutto a 0):")
        print(f"  V(buffered_delay_out) DC: {v_buffered_delay_out_dc:.3f} V")
        print(f"  V(feedback_loop_node) DC: {v_feedback_loop_node_dc:.3f} V (dovrebbe essere circa 0 * factor = 0)")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
