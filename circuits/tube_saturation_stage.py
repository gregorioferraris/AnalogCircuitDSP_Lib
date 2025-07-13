# circuits/tube_saturation_stage.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.triode import Triode # Assicurati che la tua classe Triode sia qui
from circuit_solver.circuit import Circuit

class TubeSaturationStage(Circuit):
    def __init__(self, name="Tube_Saturation_Stage"):
        super().__init__(name)

        # --- Parametri del circuito ---
        VCC = 250.0  # Tensione di alimentazione alta per valvole
        RA_val = 100e3 # Resistenza di Anodo (100k Ohm)
        RK_val = 1.5e3 # Resistenza di Catodo (1.5k Ohm) - per auto-polarizzazione
        RG_val = 1.0e6 # Resistenza di Griglia (1M Ohm) - per isolare l'ingresso DC

        C_in_val = 0.1e-6 # Condensatore di accoppiamento in ingresso (0.1uF)
        C_out_val = 0.1e-6 # Condensatore di accoppiamento in uscita (0.1uF)
        CK_val = 22e-6 # Condensatore di bypass Catodo (22uF)

        # --- Parametri della Trioide (es. 12AX7 approssimata) ---
        # Useremo i parametri tipici che hai nella tua classe Triode
        triode_params = {
            "mu": 100.0,    # Fattore di amplificazione
            "rp": 62.5e3,   # Resistenza di placca (in Ohm)
            "k": 1000.0     # Parametro di non-linearità (più alto = più non lineare)
        }

        # --- Nodi del circuito ---
        self.add_node('in')
        self.add_node('out')
        self.add_node('VCC_node') # Usiamo 'VCC_node' per evitare conflitti con il nome della sorgente
        self.add_node('grid')
        self.add_node('anode')
        self.add_node('cathode')

        # --- Componenti ---

        # 1. Alimentazione DC
        self.add_component(VoltageSource(name="VCC_Power", nodes={'pos': 'VCC_node', 'neg': 'gnd'}, voltage=VCC))

        # 2. Resistenza di Anodo (RA)
        self.add_component(Resistor(name="RA", node1='VCC_node', node2='anode', resistance=RA_val))

        # 3. Resistenza di Catodo (RK)
        self.add_component(Resistor(name="RK", node1='cathode', node2='gnd', resistance=RK_val))

        # 4. Resistenza di Griglia (RG)
        self.add_component(Resistor(name="RG", node1='grid', node2='gnd', resistance=RG_val))

        # 5. Condensatore di accoppiamento in ingresso (C_in)
        self.add_component(Capacitor(name="C_in", node1='in', node2='grid', capacitance=C_in_val))

        # 6. Condensatore di accoppiamento in uscita (C_out)
        self.add_component(Capacitor(name="C_out", node1='anode', node2='out', capacitance=C_out_val))

        # 7. Condensatore di bypass sul Catodo (CK)
        self.add_component(Capacitor(name="CK", node1='cathode', node2='gnd', capacitance=CK_val))

        # 8. Trioide
        self.add_component(Triode(name="V1",
                                  nodes={'anode': 'anode', 'grid': 'grid', 'cathode': 'cathode'},
                                  **triode_params))

        # --- Input (sorgente di segnale AC per analisi transitoria) ---
        # Per l'analisi DC, V_in_signal è 0.
        # Per vedere la saturazione, dovrai implementare una simulazione transitoria
        # e applicare un segnale d'ingresso sinusoidale con ampiezza crescente.
        self.add_component(VoltageSource(name="V_in_signal", nodes={'pos': 'in', 'neg': 'gnd'},
                                         voltage=0.0))

if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    saturation_circuit = TubeSaturationStage()
    
    print("Nodi del circuito:", saturation_circuit.get_node_names())
    print("Componenti del circuito:")
    for comp in saturation_circuit.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(saturation_circuit)

    print("\nInizio simulazione DC per trovare il punto di riposo...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in saturation_circuit.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        v1_triode = None
        for comp in saturation_circuit.nonlinear_components:
            if isinstance(comp, Triode) and comp.name == "V1":
                v1_triode = comp
                break
        
        if v1_triode:
            V_a = node_voltages[solver._get_node_index(v1_triode.nodes['anode'])]
            V_g = node_voltages[solver._get_node_index(v1_triode.nodes['grid'])]
            V_k = node_voltages[solver._get_node_index(v1_triode.nodes['cathode'])]
            
            V_gk_q = V_g - V_k
            V_ak_q = V_a - V_k
            
            Ia_q = v1_triode.calculate_plate_current(V_gk_q, V_ak_q)
            
            print(f"\n--- Punto di Riposo della Trioide (V1) ---")
            print(f"  V_GK (V1) = {V_gk_q:.3f} V")
            print(f"  V_AK (V1) = {V_ak_q:.3f} V")
            print(f"  I_A (V1) = {Ia_q*1e3:.3f} mA")
            
            # Verificare che la valvola sia in conduzione e non in cutoff profondo
            if Ia_q > 0.05e-3: # Almeno 0.05mA per considerarla attiva
                print(f"\nLa Trioide (V1) è correttamente polarizzata.")
            else:
                print(f"\nATTENZIONE: La Trioide (V1) potrebbe essere in cutoff o quasi spenta. Controlla la polarizzazione.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
