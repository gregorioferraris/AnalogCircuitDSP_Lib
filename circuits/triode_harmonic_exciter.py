# circuits/triode_harmonic_exciter.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.triode import Triode # Assicurati che la tua classe Triode sia importabile
from circuit_solver.circuit import Circuit

class TriodeHarmonicExciter(Circuit):
    def __init__(self, name="Triode_Harmonic_Exciter"):
        super().__init__(name)

        # --- Parametri del circuito ---
        V_Bplus_val = 200.0  # Volts - Alimentazione anodica
        RG_val = 1.0e6       # 1 MOhm - Resistenza di griglia
        RK_val = 1.0e3       # 1 kOhm - Resistenza di catodo
        RPlate_val = 100e3   # 100 kOhm - Resistenza anodica

        C_in_val = 0.022e-6  # 0.022 uF - Condensatore di accoppiamento in ingresso
        C_out_val = 0.022e-6 # 0.022 uF - Condensatore di accoppiamento in uscita
        CK_val = 22e-6       # 22 uF - Condensatore di bypass del catodo

        # Parametri del Triodo (es. 12AX7, un triodo ad alto guadagno, comunemente usato)
        triode_params = {
            "mu": 100.0,       # Fattore di amplificazione (ad esempio per 12AX7)
            "rp": 62.5e3,      # Resistenza interna della placca (Ohm, per 12AX7)
            "gm": 0.0016      # Transconduttanza (Siemens, per 12AX7)
        }

        # --- Nodi del circuito ---
        self.add_node('in')
        self.add_node('out')
        self.add_node('V_Bplus_node')
        self.add_node('grid')        # Griglia
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
        self.add_component(Resistor(name="RG", node1='in', node2='grid', resistance=RG_val))

        # 5. Resistenza di catodo per auto-polarizzazione
        self.add_component(Resistor(name="R_K", node1='cathode', node2='gnd', resistance=RK_val))

        # 6. Condensatore di bypass del catodo
        self.add_component(Capacitor(name="C_K", node1='cathode', node2='gnd', capacitance=CK_val))

        # 7. Resistenza anodica (Plate Resistor)
        self.add_component(Resistor(name="R_Plate", node1='V_Bplus_node', node2='plate', resistance=RPlate_val))

        # 8. Triodo (la tua classe Triode)
        self.add_component(Triode(name="V1",
                                  nodes={'plate': 'plate', 'grid': 'grid', 'cathode': 'cathode'},
                                  **triode_params))

        # 9. Condensatore di accoppiamento in uscita
        self.add_component(Capacitor(name="C_out", node1='plate', node2='out', capacitance=C_out_val))
        
        # 10. Resistenza di carico in uscita
        self.add_component(Resistor(name="R_Load", node1='out', node2='gnd', resistance=100e3)) # 100k Ohm di carico

# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    triode_exciter = TriodeHarmonicExciter()
    
    print(f"Circuito: {triode_exciter.name}")
    print("Nodi:", triode_exciter.get_node_names())
    print("Componenti:")
    for comp in triode_exciter.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(triode_exciter)

    print("\nSimulazione DC del Triode Harmonic Exciter...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in triode_exciter.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Verifica il punto di riposo del Triodo
        v1 = next((c for c in triode_exciter.nonlinear_components if isinstance(c, Triode) and c.name == "V1"), None)
        if v1:
            V_plate_q = node_voltages[solver._get_node_index(v1.nodes['plate'])]
            V_grid_q = node_voltages[solver._get_node_index(v1.nodes['grid'])]
            V_cathode_q = node_voltages[solver._get_node_index(v1.nodes['cathode'])]
            
            Vgk_q = V_grid_q - V_cathode_q
            Vpk_q = V_plate_q - V_cathode_q
            Ip_q = v1.calculate_plate_current(Vgk_q, Vpk_q)

            print(f"\n--- Punto di Riposo del Triodo (V1) ---")
            print(f"  V_GK (V1) = {Vgk_q:.3f} V")
            print(f"  V_PK (V1) = {Vpk_q:.3f} V")
            print(f"  I_P (V1) = {Ip_q*1e3:.3f} mA")
            
            # Condizioni per il funzionamento in regione attiva
            if Vgk_q < 0 and Vpk_q > 50 and Ip_q > 0.01e-3: # Condizione semplificata di attiva
                print(f"  V1 è in regione Attiva (per l'amplificazione).")
            elif Vgk_q >= 0:
                print(f"  V1 è in conduzione di griglia (potrebbe non essere l'ideale per audio pulito).")
            else:
                print(f"  V1 è in Cut-off (spento) o quasi.")

            if Ip_q > 0.5e-3 and Ip_q < 10e-3: # Esempio: corrente tra 0.5mA e 10mA per un piccolo segnale
                print(f"Il Triodo (V1) è correttamente polarizzato per operare in Classe A.")
            else:
                print(f"ATTENZIONE: Il Triodo (V1) potrebbe essere fuori dalla polarizzazione desiderata.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
