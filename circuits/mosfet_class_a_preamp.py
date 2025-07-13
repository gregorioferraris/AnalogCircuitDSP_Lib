# circuits/mosfet_class_a_preamp.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.mosfet import MOSFET # Assicurati che la tua classe MOSFET sia qui
from circuit_solver.circuit import Circuit

class MOSFETClassAPreamp(Circuit):
    def __init__(self, name="MOSFET_Class_A_Preamp"):
        super().__init__(name)

        # --- Parametri del circuito ---
        VCC = 9.0  # Tensione di alimentazione
        R1_val = 1.0e6 # Resistenza divisore di tensione (Gate 1) - molto alta per MOSFET
        R2_val = 330e3 # Resistenza divisore di tensione (Gate 2) - molto alta
        RD_val = 2.2e3 # Resistenza di Drain
        RS_val = 470   # Resistenza di Source

        C_in_val = 10e-6 # Condensatore di accoppiamento in ingresso (10uF)
        C_out_val = 10e-6 # Condensatore di accoppiamento in uscita (10uF)
        CS_val = 100e-6 # Condensatore di bypass Source (100uF)

        # --- Parametri del MOSFET (N-channel Enhancement-mode generico) ---
        # Assicurati che la tua classe MOSFET abbia questi parametri o simili.
        # k è il parametro di transconduttanza (mA/V^2), Vth è la tensione di soglia (V)
        mosfet_params = {
            "k": 0.003,     # Transconductance parameter (A/V^2), e.g., 3mA/V^2
            "Vth": 1.5,     # Threshold Voltage (V)
            "lambda_val": 0.02 # Channel-length modulation parameter
        }

        # --- Nodi del circuito ---
        self.add_node('in')      # Ingresso audio
        self.add_node('out')     # Uscita audio
        self.add_node('VCC')     # Nodo di alimentazione
        self.add_node('gate')    # Gate del MOSFET
        self.add_node('drain')   # Drain del MOSFET
        self.add_node('source')  # Source del MOSFET

        # --- Componenti ---

        # 1. Alimentazione DC
        self.add_component(VoltageSource(name="VCC_Source", nodes={'pos': 'VCC', 'neg': 'gnd'}, voltage=VCC))

        # 2. Divisore di tensione per il Gate (R1, R2)
        self.add_component(Resistor(name="R1", node1='VCC', node2='gate', resistance=R1_val))
        self.add_component(Resistor(name="R2", node1='gate', node2='gnd', resistance=R2_val))

        # 3. Resistenza di Drain (RD)
        self.add_component(Resistor(name="RD", node1='VCC', node2='drain', resistance=RD_val))

        # 4. Resistenza di Source (RS)
        self.add_component(Resistor(name="RS", node1='source', node2='gnd', resistance=RS_val))

        # 5. Condensatore di accoppiamento in ingresso (C_in)
        self.add_component(Capacitor(name="C_in", node1='in', node2='gate', capacitance=C_in_val))

        # 6. Condensatore di accoppiamento in uscita (C_out)
        self.add_component(Capacitor(name="C_out", node1='drain', node2='out', capacitance=C_out_val))

        # 7. Condensatore di bypass sul Source (CS)
        self.add_component(Capacitor(name="CS", node1='source', node2='gnd', capacitance=CS_val))

        # 8. Transistor MOSFET N-channel Enhancement-mode
        # Assicurati che la tua classe MOSFET supporti l'inizializzazione con k, Vth, lambda_val.
        # E che abbia una funzione calculate_drain_current(Vgs, Vds) e calculate_jacobian_elements.
        self.add_component(MOSFET(name="M1",
                                nodes={'drain': 'drain', 'gate': 'gate', 'source': 'source'},
                                **mosfet_params))

        # --- Input (sorgente di segnale AC) ---
        self.add_component(VoltageSource(name="V_in_signal", nodes={'pos': 'in', 'neg': 'gnd'},
                                         voltage=0.0)) # Valor 0 per DC analysis.

if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    preamp_circuit = MOSFETClassAPreamp()
    
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
        
        m1_mosfet = None
        for comp in preamp_circuit.nonlinear_components:
            if isinstance(comp, MOSFET) and comp.name == "M1":
                m1_mosfet = comp
                break
        
        if m1_mosfet:
            V_d = node_voltages[solver._get_node_index(m1_mosfet.nodes['drain'])]
            V_g = node_voltages[solver._get_node_index(m1_mosfet.nodes['gate'])]
            V_s = node_voltages[solver._get_node_index(m1_mosfet.nodes['source'])]
            
            V_gs_q = V_g - V_s
            V_ds_q = V_d - V_s
            
            Id_q = m1_mosfet.calculate_drain_current(V_gs_q, V_ds_q)
            
            print(f"\n--- Punto di Riposo del MOSFET (M1) ---")
            print(f"  V_GS (M1) = {V_gs_q:.3f} V")
            print(f"  V_DS (M1) = {V_ds_q:.3f} V")
            print(f"  I_D (M1) = {Id_q*1e3:.3f} mA")
            
            # Controllo approssimato per la regione di saturazione del MOSFET Enhancement-mode
            if Id_q > 0.0 and V_gs_q > m1_mosfet.Vth and V_ds_q >= (V_gs_q - m1_mosfet.Vth):
                print(f"\nIl MOSFET (M1) è polarizzato in Classe A nella regione di saturazione.")
            else:
                print(f"\nATTENZIONE: Il MOSFET (M1) potrebbe non essere polarizzato correttamente in Classe A (Id = 0, Vgs <= Vth, o in triodo).")
                print("Controlla i valori delle resistenze di polarizzazione e la tensione Vcc.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
