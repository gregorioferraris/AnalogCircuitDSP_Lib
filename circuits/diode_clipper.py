# circuits/diode_clipper.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.diode import Diode # Assicurati che la tua classe Diode sia qui
from circuit_solver.circuit import Circuit

class DiodeClipper(Circuit):
    def __init__(self, name="Diode_Clipper"):
        super().__init__(name)

        # --- Parametri del circuito ---
        # Alimentazione può essere singola o duale; qui usiamo singola per semplicità.
        # Per un clipping simmetrico attorno a 0, servirebbe un bias o alimentazione duale.
        # Qui un clipping che avviene sopra 0.7V e sotto -0.7V (se collegati a gnd).
        VCC = 9.0 # Alimentazione per un eventuale buffer o solo per testare i diodi.
        # Per la pura rete di clipping, potremmo anche non aver bisogno di una VCC esplicita,
        # dato che i diodi si basano solo sulla tensione differenziale.

        R_in_val = 10e3 # Resistenza di ingresso (10k Ohm)
        R_out_val = 10e3 # Resistenza in serie ai diodi (10k Ohm)
        
        C_in_val = 10e-6 # Condensatore di accoppiamento in ingresso
        C_out_val = 10e-6 # Condensatore di accoppiamento in uscita

        # --- Parametri del Diodo (es. 1N4148 approssimato) ---
        diode_params = {
            "Is": 1.0e-14, # Corrente di saturazione inversa (A)
            "Vt": 0.02585, # Tensione termica a 300K (V)
            "N": 1.0       # Fattore di idealità
        }

        # --- Nodi del circuito ---
        self.add_node('in')
        self.add_node('clip_node') # Nodo dove avviene il clipping
        self.add_node('out')
        self.add_node('VCC_node') # Per l'alimentazione, anche se non usata direttamente dai diodi

        # --- Componenti ---

        # 1. Alimentazione (se necessaria per bias o altri stadi)
        self.add_component(VoltageSource(name="VCC_Power", nodes={'pos': 'VCC_node', 'neg': 'gnd'}, voltage=VCC))

        # 2. Resistenza di ingresso (per limitare la corrente verso i diodi)
        self.add_component(Resistor(name="R_in", node1='in', node2='clip_node', resistance=R_in_val))

        # 3. Diodo 1: dall'anodo al gnd, catodo al clip_node (conduce per V_clip_node negativa)
        # Questo creerà clipping per la parte negativa del segnale.
        self.add_component(Diode(name="D1", nodes={'anode': 'gnd', 'cathode': 'clip_node'}, **diode_params))

        # 4. Diodo 2: dall'anodo al clip_node, catodo al gnd (conduce per V_clip_node positiva)
        # Questo creerà clipping per la parte positiva del segnale.
        self.add_component(Diode(name="D2", nodes={'anode': 'clip_node', 'cathode': 'gnd'}, **diode_params))
        
        # 5. Resistenza in uscita (o carico)
        self.add_component(Resistor(name="R_out", node1='clip_node', node2='out', resistance=R_out_val))

        # 6. Condensatore di accoppiamento in ingresso
        self.add_component(Capacitor(name="C_in", node1='V_in_signal_pos', node2='in', capacitance=C_in_val))

        # 7. Condensatore di accoppiamento in uscita
        self.add_component(Capacitor(name="C_out", node1='out', node2='gnd', capacitance=C_out_val))

        # --- Input (sorgente di segnale AC per analisi transitoria) ---
        # Per vedere il clipping, devi applicare un segnale sinusoidale in AC/transitorio
        # e osservare la forma d'onda in 'out'. Per DC, è 0V.
        self.add_component(VoltageSource(name="V_in_signal", nodes={'pos': 'V_in_signal_pos', 'neg': 'gnd'},
                                         voltage=0.0))
        # Aggiungo un nodo V_in_signal_pos per l'ingresso AC per distinguerlo dal nodo 'in' dopo il C_in
        self.add_node('V_in_signal_pos')

if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    clipper_circuit = DiodeClipper()
    
    print("Nodi del circuito:", clipper_circuit.get_node_names())
    print("Componenti del circuito:")
    for comp in clipper_circuit.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(clipper_circuit)

    print("\nInizio simulazione DC per il distorsore a diodi...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in clipper_circuit.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Per il clipper a diodi, il Q-point non è così critico come per gli amplificatori,
        # dato che i diodi sono passivi e si attivano solo con il segnale.
        # Ma possiamo controllare lo stato dei diodi.
        d1_diode = None
        d2_diode = None
        for comp in clipper_circuit.nonlinear_components:
            if isinstance(comp, Diode):
                if comp.name == "D1": d1_diode = comp
                if comp.name == "D2": d2_diode = comp
        
        if d1_diode and d2_diode:
            V_clip_node_q = node_voltages[solver._get_node_index('clip_node')]
            
            Id1_q = d1_diode.calculate_diode_current(0.0 - V_clip_node_q) # Anodo=gnd, Catodo=clip_node
            Id2_q = d2_diode.calculate_diode_current(V_clip_node_q - 0.0) # Anodo=clip_node, Catodo=gnd
            
            print(f"\n--- Correnti dei Diodi al Punto di Riposo ---")
            print(f"  I_D1 = {Id1_q*1e6:.3f} uA")
            print(f"  I_D2 = {Id2_q*1e6:.3f} uA")
            print(f"  V_clip_node (DC) = {V_clip_node_q:.3f} V")
            
            # In DC senza segnale, le correnti dovrebbero essere molto basse (vicino allo zero).
            if abs(Id1_q) < 1e-9 and abs(Id2_q) < 1e-9:
                print("\nI diodi sono spenti al punto di riposo DC, come previsto per un clipper passivo.")
            else:
                print("\nATTENZIONE: I diodi potrebbero avere una corrente significativa al punto di riposo DC.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
