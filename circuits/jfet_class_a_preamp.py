# circuits/jfet_class_a_preamp.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.jfet import JFET # Assicurati che la tua classe JFET sia qui
from circuit_solver.circuit import Circuit

class JFETClassAPreamp(Circuit):
    def __init__(self, name="JFET_Class_A_Preamp"):
        super().__init__(name)

        # --- Parametri del circuito ---
        VCC = 9.0  # Tensione di alimentazione
        RD_val = 2.2e3 # Resistenza di Drain
        RS_val = 1.0e3 # Resistenza di Source (per auto-polarizzazione)

        C_in_val = 10e-6 # Condensatore di accoppiamento in ingresso (10uF)
        C_out_val = 10e-6 # Condensatore di accoppiamento in uscita (10uF)
        CS_val = 100e-6 # Condensatore di bypass Source (100uF)

        # --- Parametri del JFET (N-channel generico, es. J201 approssimato) ---
        jfet_params = {
            "Idss": 0.003,  # 3 mA (Corrente di Drain a Vgs=0)
            "Vp": -0.8,     # -0.8 V (Tensione di Pinch-Off)
            "lambda_val": 0.01 # Piccolo effetto di modulazione di canale
        }

        # --- Nodi del circuito ---
        self.add_node('in')      # Ingresso audio
        self.add_node('out')     # Uscita audio
        self.add_node('VCC')     # Nodo di alimentazione
        self.add_node('gate')    # Gate del JFET
        self.add_node('drain')   # Drain del JFET
        self.add_node('source')  # Source del JFET

        # --- Componenti ---

        # 1. Alimentazione DC
        self.add_component(VoltageSource(name="VCC_Source", nodes={'pos': 'VCC', 'neg': 'gnd'}, voltage=VCC))

        # 2. Resistenza di Drain (RD)
        self.add_component(Resistor(name="RD", node1='VCC', node2='drain', resistance=RD_val))

        # 3. Resistenza di Source (RS) - per auto-polarizzazione
        self.add_component(Resistor(name="RS", node1='source', node2='gnd', resistance=RS_val))
        
        # 4. Resistenza di Gate (RG) - per mantenere il gate a gnd DC, ma ad alta impedenza per il segnale
        # Spesso RG è molto grande (es. 1M Ohm) per non caricare la sorgente.
        # Collega il gate a gnd (DC) ma lascia il gate "fluttuare" per il segnale AC.
        self.add_component(Resistor(name="RG", node1='gate', node2='gnd', resistance=1.0e6))

        # 5. Condensatore di accoppiamento in ingresso (C_in)
        self.add_component(Capacitor(name="C_in", node1='in', node2='gate', capacitance=C_in_val))

        # 6. Condensatore di accoppiamento in uscita (C_out)
        self.add_component(Capacitor(name="C_out", node1='drain', node2='out', capacitance=C_out_val))

        # 7. Condensatore di bypass sul Source (CS)
        self.add_component(Capacitor(name="CS", node1='source', node2='gnd', capacitance=CS_val))

        # 8. Transistor JFET N-channel
        self.add_component(JFET(name="J1",
                               nodes={'drain': 'drain', 'gate': 'gate', 'source': 'source'},
                               **jfet_params))

        # --- Input (sorgente di segnale AC per la simulazione transitoria) ---
        self.add_component(VoltageSource(name="V_in_signal", nodes={'pos': 'in', 'neg': 'gnd'},
                                         voltage=0.0)) # Valor 0 per DC analysis.

if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    preamp_circuit = JFETClassAPreamp()
    
    print("Nodi del circuito:", preamp_circuit.get_node_names())
    print("Componenti del circuito:")
    for comp in preamp_circuit.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(preamp_circuit)

    print("\nInizio simulazione DC per trovare il punto di riposo (Q-point)...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in preamp_circuit.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        j1_jfet = None
        for comp in preamp_circuit.nonlinear_components:
            if isinstance(comp, JFET) and comp.name == "J1":
                j1_jfet = comp
                break
        
        if j1_jfet:
            V_d = node_voltages[solver._get_node_index(j1_jfet.nodes['drain'])]
            V_g = node_voltages[solver._get_node_index(j1_jfet.nodes['gate'])]
            V_s = node_voltages[solver._get_node_index(j1_jfet.nodes['source'])]
            
            V_gs_q = V_g - V_s
            V_ds_q = V_d - V_s
            
            Id_q = j1_jfet.calculate_drain_current(V_gs_q, V_ds_q)
            
            print(f"\n--- Punto di Riposo del JFET (J1) ---")
            print(f"  V_GS (J1) = {V_gs_q:.3f} V")
            print(f"  V_DS (J1) = {V_ds_q:.3f} V")
            print(f"  I_D (J1) = {Id_q*1e3:.3f} mA")
            
            if Id_q > 0.0 and V_ds_q >= (V_gs_q - j1_jfet.Vp): # Controllo approssimato per la regione di saturazione
                print(f"\nIl JFET (J1) è polarizzato in Classe A nella regione di saturazione.")
            else:
                print(f"\nATTENZIONE: Il JFET (J1) potrebbe non essere polarizzato correttamente in Classe A (Id = 0 o in triodo).")
                print("Controlla i valori delle resistenze di polarizzazione e la tensione Vcc.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
