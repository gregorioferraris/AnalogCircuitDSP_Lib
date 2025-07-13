# circuits/analog_buffer.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.mosfet import MOSFET
from circuit_solver.circuit import Circuit

class AnalogBuffer(Circuit):
    def __init__(self, name="Analog_Buffer"):
        super().__init__(name)

        # --- Parametri del circuito ---
        VCC_val = 9.0      # Alimentazione
        RG_val = 1.0e6     # Resistenza di Griglia (per bias, molto alta)
        RS_val = 10e3      # Resistenza di Source (carico del follower)

        C_in_val = 1.0e-6  # Condensatore di accoppiamento in ingresso
        C_out_val = 1.0e-6 # Condensatore di accoppiamento in uscita

        # Parametri del MOSFET (ad esempio, un JFET o MOSFET a canale N come 2N7000 o BS170)
        # Importante: Vth dovrebbe essere inferiore a V(gate)-V(source) in DC per operare.
        mosfet_params = {
            "k": 0.003,      # Transconduttanza (mA/V^2)
            "Vth": 1.5,      # Tensione di soglia (V)
            "lambda_val": 0.02 # Modulazione della lunghezza del canale
        }

        # --- Nodi del circuito ---
        self.add_node('in')           # Ingresso AC del buffer
        self.add_node('out')          # Uscita AC del buffer
        self.add_node('VCC_node')     # Nodo di alimentazione
        self.add_node('gate')         # Gate del MOSFET
        self.add_node('source')       # Source del MOSFET
        self.add_node('drain')        # Drain del MOSFET (collegato a VCC per Source Follower)
        self.add_node('V_in_AC')      # Nodo per la sorgente di segnale AC

        # --- Componenti ---

        # Alimentazione DC
        self.add_component(VoltageSource(name="VCC_Power", nodes={'pos': 'VCC_node', 'neg': 'gnd'}, voltage=VCC_val))

        # Sorgente di segnale d'ingresso AC (per analisi transitoria)
        self.add_component(VoltageSource(name="V_signal", nodes={'pos': 'V_in_AC', 'neg': 'gnd'}, voltage=0.0))

        # Condensatore di accoppiamento in ingresso
        self.add_component(Capacitor(name="C_in", node1='V_in_AC', node2='in', capacitance=C_in_val))

        # Resistenza di Griglia (bias a gnd per DC, alta impedenza)
        self.add_component(Resistor(name="RG", node1='in', node2='gate', resistance=RG_val))

        # MOSFET: Drain a VCC, Source a RS, Gate all'ingresso (dopo C_in e RG)
        self.add_component(MOSFET(name="M1_buffer",
                                  nodes={'drain': 'VCC_node', 'gate': 'gate', 'source': 'source'},
                                  **mosfet_params))

        # Resistenza di Source (carico)
        self.add_component(Resistor(name="RS", node1='source', node2='gnd', resistance=RS_val))

        # Condensatore di accoppiamento in uscita
        self.add_component(Capacitor(name="C_out", node1='source', node2='out', capacitance=C_out_val))
        
        # Resistenza di carico per l'uscita del buffer
        self.add_component(Resistor(name="R_load_out", node1='out', node2='gnd', resistance=10e3))

# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    buffer_circuit = AnalogBuffer()
    
    print(f"Circuito: {buffer_circuit.name}")
    print("Nodi:", buffer_circuit.get_node_names())
    print("Componenti:")
    for comp in buffer_circuit.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(buffer_circuit)

    print("\nSimulazione DC dell'Analog Buffer...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in buffer_circuit.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        m1 = next((c for c in buffer_circuit.nonlinear_components if isinstance(c, MOSFET) and c.name == "M1_buffer"), None)
        if m1:
            V_gate_q = node_voltages[solver._get_node_index(m1.nodes['gate'])]
            V_source_q = node_voltages[solver._get_node_index(m1.nodes['source'])]
            Vgs_q = V_gate_q - V_source_q
            print(f"\n  V_GS (M1_buffer) = {Vgs_q:.3f} V")
            if Vgs_q > m1.Vth:
                print("  Il MOSFET è correttamente in conduzione.")
                # Per un buffer, ci aspettiamo V(out)DC molto vicino a 0V se C_out c'è
                print(f"  V(out) DC: {node_voltages[solver._get_node_index('out')]:.3f} V")
            else:
                print("ATTENZIONE: Il MOSFET potrebbe non essere polarizzato correttamente (cutoff).")
    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
