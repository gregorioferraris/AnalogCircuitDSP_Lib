# circuits/bjt_class_a_preamp.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.bjt import BJT # Assicurati che la tua classe BJT sia qui
from circuit_solver.circuit import Circuit

class BJTClassAPreamp(Circuit):
    def __init__(self, name="BJT_Class_A_Preamp"):
        super().__init__(name)

        # --- Parametri del circuito ---
        VCC = 9.0  # Tensione di alimentazione
        R1_val = 100e3 # Resistenza divisore di tensione (Base 1)
        R2_val = 22e3  # Resistenza divisore di tensione (Base 2)
        RC_val = 3.3e3 # Resistenza di Collettore
        RE_val = 1.0e3 # Resistenza di Emettitore

        C_in_val = 10e-6 # Condensatore di accoppiamento in ingresso (10uF)
        C_out_val = 10e-6 # Condensatore di accoppiamento in uscita (10uF)
        CE_val = 100e-6 # Condensatore di bypass emettitore (100uF) - per un guadagno AC più alto

        # --- Parametri del BJT (NPN generico, es. 2N3904 approssimato) ---
        # Useremo i parametri tipici che hai nella tua classe BJT
        # Se la tua classe BJT richiede un 'name' e 'nodes' nel costruttore, adegua.
        bjt_params = {
            "Is": 1.0e-14,    # A (Corrente di saturazione inversa)
            "Vt": 0.02585,    # V (Tensione termica a 300K)
            "Beta_F": 200.0,  # hFE (Guadagno di corrente in polarizzazione diretta)
            "Va": 70.0        # V (Tensione di Early)
        }

        # --- Nodi del circuito ---
        # GROUND è implicitamente il nodo 0 o il nodo di riferimento del solutore
        # Noi useremo 'gnd' come nome convenzionale
        self.add_node('in')      # Ingresso audio
        self.add_node('out')     # Uscita audio
        self.add_node('VCC')     # Nodo di alimentazione
        self.add_node('base')    # Base del BJT
        self.add_node('collector') # Collettore del BJT
        self.add_node('emitter') # Emettitore del BJT

        # --- Componenti ---

        # 1. Alimentazione DC
        self.add_component(VoltageSource(name="VCC_Source", nodes={'pos': 'VCC', 'neg': 'gnd'}, voltage=VCC))

        # 2. Divisore di tensione per la base (R1, R2)
        self.add_component(Resistor(name="R1", node1='VCC', node2='base', resistance=R1_val))
        self.add_component(Resistor(name="R2", node1='base', node2='gnd', resistance=R2_val))

        # 3. Resistenza di Collettore (RC)
        self.add_component(Resistor(name="RC", node1='VCC', node2='collector', resistance=RC_val))

        # 4. Resistenza di Emettitore (RE)
        self.add_component(Resistor(name="RE", node1='emitter', node2='gnd', resistance=RE_val))

        # 5. Condensatore di accoppiamento in ingresso (C_in)
        self.add_component(Capacitor(name="C_in", node1='in', node2='base', capacitance=C_in_val))

        # 6. Condensatore di accoppiamento in uscita (C_out)
        self.add_component(Capacitor(name="C_out", node1='collector', node2='out', capacitance=C_out_val))

        # 7. Condensatore di bypass sull'emettitore (CE)
        self.add_component(Capacitor(name="CE", node1='emitter', node2='gnd', capacitance=CE_val))

        # 8. Transistor BJT NPN
        # La classe BJT dovrebbe prendere i nomi dei nodi così:
        self.add_component(BJT(name="Q1",
                               nodes={'collector': 'collector', 'base': 'base', 'emitter': 'emitter'},
                               **bjt_params))

        # --- Input (sorgente di segnale AC per la simulazione transitoria) ---
        # Durante la simulazione DC, V_in sarà 0.
        # Durante la simulazione AC (se la implementerai), userai questa sorgente.
        # Per la DC, la sorgente AC ha valore 0.
        self.add_component(VoltageSource(name="V_in_signal", nodes={'pos': 'in', 'neg': 'gnd'},
                                         voltage=0.0)) # Valor 0 per DC analysis. Per AC/Transient sarà una sinusoide.

if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver # Assicurati che il percorso sia corretto

    # Crea un'istanza del preamplificatore
    preamp_circuit = BJTClassAPreamp()
    
    # Stampa i nodi e i componenti per verifica
    print("Nodi del circuito:", preamp_circuit.get_node_names())
    print("Componenti del circuito:")
    for comp in preamp_circuit.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    # Crea il solutore MNA
    solver = MnaSolver(preamp_circuit)

    print("\nInizio simulazione DC per trovare il punto di riposo (Q-point)...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        # Stampa le tensioni ai nodi in modo leggibile
        print("Tensioni ai Nodi:")
        for node_name, idx in preamp_circuit.get_node_map().items():
            if idx != -1: # Ignora il nodo Ground se ha indice -1
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Estrai le correnti di interesse dal BJT al Q-point
        # Avrai bisogno di accedere al BJT specifico dal circuito
        q1_bjt = None
        for comp in preamp_circuit.nonlinear_components:
            if isinstance(comp, BJT) and comp.name == "Q1":
                q1_bjt = comp
                break
        
        if q1_bjt:
            # Recupera le tensioni DC ai nodi del BJT
            V_c = node_voltages[solver._get_node_index(q1_bjt.nodes['collector'])]
            V_b = node_voltages[solver._get_node_index(q1_bjt.nodes['base'])]
            V_e = node_voltages[solver._get_node_index(q1_bjt.nodes['emitter'])]
            
            V_be_q = V_b - V_e
            V_ce_q = V_c - V_e
            
            Ic_q = q1_bjt.calculate_collector_current(V_be_q, V_ce_q)
            Ib_q = q1_bjt.calculate_base_current(Ic_q)
            
            print(f"\n--- Punto di Riposo del BJT (Q1) ---")
            print(f"  V_BE (Q1) = {V_be_q:.3f} V")
            print(f"  V_CE (Q1) = {V_ce_q:.3f} V")
            print(f"  I_C (Q1) = {Ic_q*1e3:.3f} mA")
            print(f"  I_B (Q1) = {Ib_q*1e6:.3f} uA")
            
            # Valutazione rapida della Classe A
            if Ic_q > 0.0:
                print(f"\nIl BJT (Q1) è polarizzato in Classe A (Ic > 0).")
            else:
                print(f"\nATTENZIONE: Il BJT (Q1) potrebbe non essere polarizzato correttamente in Classe A (Ic = 0).")
                print("Controlla i valori delle resistenze di polarizzazione e la tensione Vcc.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
