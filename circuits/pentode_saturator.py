# circuits/pentode_saturator.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.pentode import Pentode # Assicurati che la tua classe Pentode sia importabile
from circuit_solver.circuit import Circuit

class PentodeSaturator(Circuit):
    def __init__(self, name="Pentode_Saturator"):
        super().__init__(name)

        # --- Parametri del circuito ---
        V_Bplus_val = 250.0  # Volts - Alimentazione anodica (tipica per valvole audio)
        RG1_val = 1.0e6      # 1 MOhm - Resistenza di griglia (bias di ingresso)
        RK_val = 1.5e3       # 1.5 kOhm - Resistenza di catodo (per auto-polarizzazione)
        RScreen_val = 100e3  # 100 kOhm - Resistenza per la griglia schermo (limita corrente, imposta tensione)
        RPlate_val = 47e3    # 47 kOhm - Resistenza anodica (carico)

        C_in_val = 0.1e-6    # 0.1 uF - Condensatore di accoppiamento in ingresso
        C_out_val = 0.1e-6   # 0.1 uF - Condensatore di accoppiamento in uscita
        CK_val = 22e-6       # 22 uF - Condensatore di bypass del catodo
        CScreen_val = 0.1e-6 # 0.1 uF - Condensatore di bypass della griglia schermo

        # Parametri del Pentodo (es. EL84 approssimato per preamp-like usage)
        # Nota: I parametri per un pentodo sono più complessi e variano molto.
        # Questi sono valori di esempio.
        pentode_params = {
            "mu": 100.0,       # Fattore di amplificazione (valore alto per pentodi)
            "rp": 50e3,        # Resistenza interna della placca (Ohm)
            "gm": 0.005,       # Transconduttanza (Siemens) - Questo è fondamentale per la corrente di placca
            "Vth": -2.0,       # Tensione di soglia della griglia di controllo (V)
            "k_saturation": 1.5e-3 # Parametro per la saturazione (mA/V^2) - adatta la formula
        }
        # In un modello Pentodo più robusto, avresti anche la corrente di griglia schermo.
        # Per questa simulazione semplificata, ci concentriamo sulla corrente anodica.

        # --- Nodi del circuito ---
        self.add_node('in')
        self.add_node('out')
        self.add_node('V_Bplus_node')
        self.add_node('g1_grid')     # Griglia di controllo (G1)
        self.add_node('screen_grid') # Griglia schermo (G2)
        self.add_node('plate')       # Anodo / Placca
        self.add_node('cathode')     # Catodo
        self.add_node('V_in_AC')     # Nodo per la sorgente di segnale AC

        # --- Componenti ---

        # 1. Alimentazione DC
        self.add_component(VoltageSource(name="V_Bplus", nodes={'pos': 'V_Bplus_node', 'neg': 'gnd'}, voltage=V_Bplus_val))

        # 2. Sorgente di segnale d'ingresso AC
        self.add_component(VoltageSource(name="V_signal", nodes={'pos': 'V_in_AC', 'neg': 'gnd'}, voltage=0.0))

        # 3. Condensatore di accoppiamento in ingresso
        self.add_component(Capacitor(name="C_in", node1='V_in_AC', node2='in', capacitance=C_in_val))

        # 4. Resistenza di griglia (Grid Leak Resistor)
        self.add_component(Resistor(name="RG1", node1='in', node2='g1_grid', resistance=RG1_val))

        # 5. Resistenza di catodo per auto-polarizzazione
        self.add_component(Resistor(name="R_K", node1='cathode', node2='gnd', resistance=RK_val))

        # 6. Condensatore di bypass del catodo
        self.add_component(Capacitor(name="C_K", node1='cathode', node2='gnd', capacitance=CK_val))

        # 7. Resistenza della griglia schermo (limita corrente G2 e imposta tensione G2)
        self.add_component(Resistor(name="R_Screen", node1='V_Bplus_node', node2='screen_grid', resistance=RScreen_val))

        # 8. Condensatore di bypass della griglia schermo (a massa)
        self.add_component(Capacitor(name="C_Screen", node1='screen_grid', node2='gnd', capacitance=CScreen_val))

        # 9. Resistenza anodica (Plate Resistor)
        self.add_component(Resistor(name="R_Plate", node1='V_Bplus_node', node2='plate', resistance=RPlate_val))

        # 10. Pentodo (la tua classe Pentode)
        self.add_component(Pentode(name="V1",
                                    nodes={'plate': 'plate', 'g1_grid': 'g1_grid', 'screen_grid': 'screen_grid', 'cathode': 'cathode'},
                                    **pentode_params))

        # 11. Condensatore di accoppiamento in uscita
        self.add_component(Capacitor(name="C_out", node1='plate', node2='out', capacitance=C_out_val))

        # 12. Resistenza di carico in uscita
        self.add_component(Resistor(name="R_Load", node1='out', node2='gnd', resistance=100e3)) # 100k Ohm di carico

# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    pentode_saturator = PentodeSaturator()
    
    print(f"Circuito: {pentode_saturator.name}")
    print("Nodi:", pentode_saturator.get_node_names())
    print("Componenti:")
    for comp in pentode_saturator.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(pentode_saturator)

    print("\nSimulazione DC del Pentode Saturator...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in pentode_saturator.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Verifica il punto di riposo del Pentodo
        v1 = next((c for c in pentode_saturator.nonlinear_components if isinstance(c, Pentode) and c.name == "V1"), None)
        if v1:
            V_plate_q = node_voltages[solver._get_node_index(v1.nodes['plate'])]
            V_g1_grid_q = node_voltages[solver._get_node_index(v1.nodes['g1_grid'])]
            V_screen_grid_q = node_voltages[solver._get_node_index(v1.nodes['screen_grid'])]
            V_cathode_q = node_voltages[solver._get_node_index(v1.nodes['cathode'])]
            
            Vgk_q = V_g1_grid_q - V_cathode_q
            Vpk_q = V_plate_q - V_cathode_q
            Vsg_q = V_screen_grid_q - V_cathode_q # Tensione screen-catodo
            
            Ip_q = v1.calculate_plate_current(Vgk_q, Vpk_q, Vsg_q)

            print(f"\n--- Punto di Riposo del Pentodo (V1) ---")
            print(f"  V_GK (V1) = {Vgk_q:.3f} V")
            print(f"  V_PK (V1) = {Vpk_q:.3f} V")
            print(f"  V_SG (V1) = {Vsg_q:.3f} V")
            print(f"  I_P (V1) = {Ip_q*1e3:.3f} mA")
            
            # Condizioni per il funzionamento in regione attiva (per amplificazione)
            # In saturazione, Ip si avvicina al limite massimo.
            if Vgk_q < v1.Vth and Vpk_q > 50 and Ip_q > 0.01e-3: # Condizione semplificata di attiva
                print(f"  V1 è in regione Attiva (per l'amplificazione).")
            elif Vgk_q >= v1.Vth:
                print(f"  V1 è in Cut-off (spento).")
            else:
                print(f"  V1 è in Triodo o Saturazione (Ip alta).")

            if Ip_q > 0.5e-3 and Ip_q < 20e-3: # Esempio: corrente tra 0.5mA e 20mA per un piccolo segnale/preamp
                print(f"Il Pentodo (V1) è correttamente polarizzato per operare in Classe A.")
            else:
                print(f"ATTENZIONE: Il Pentodo (V1) potrebbe essere fuori dalla polarizzazione desiderata.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
