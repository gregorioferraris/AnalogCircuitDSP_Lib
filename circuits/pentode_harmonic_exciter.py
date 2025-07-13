# circuits/pentode_harmonic_exciter.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.pentode import Pentode # Assicurati che la tua classe Pentode sia importabile
from circuit_solver.circuit import Circuit

class PentodeHarmonicExciter(Circuit):
    def __init__(self, name="Pentode_Harmonic_Exciter"):
        super().__init__(name)

        # --- Parametri del circuito ---
        # Si usano tensioni e resistenze simili a quelle di un amplificatore per chitarra/preamp
        V_Bplus_val = 300.0  # Volts - Alimentazione anodica più alta per un maggiore swing
        RG1_val = 1.0e6      # 1 MOhm
        RK_val = 1.0e3       # 1 kOhm - Potrebbe essere variato per alterare il bias e la distorsione
        RScreen_val = 50e3   # 50 kOhm - Resistenza per la griglia schermo
        RPlate_val = 22e3    # 22 kOhm - Resistenza anodica (più bassa per maggiore corrente e distorsione)

        C_in_val = 0.047e-6  # 0.047 uF
        C_out_val = 0.047e-6 # 0.047 uF
        CK_val = 22e-6       # 22 uF
        CScreen_val = 0.47e-6 # 0.47 uF - Condensatore di bypass più grande per stabilità G2

        # Parametri del Pentodo (es. EL84 o 6V6 approssimato, più votato alla distorsione)
        pentode_params = {
            "mu": 100.0,       # Fattore di amplificazione
            "rp": 25e3,        # Resistenza interna della placca (Ohm, leggermente più bassa per saturazione)
            "gm": 0.005,       # Transconduttanza (Siemens)
            "Vth": -3.0,       # Tensione di soglia della griglia di controllo (V) - potremmo renderla più negativa
            "k_saturation": 2.0e-3 # Parametro di saturazione, adattato per una maggiore non-linearità
        }

        # --- Nodi del circuito ---
        self.add_node('in')
        self.add_node('out')
        self.add_node('V_Bplus_node')
        self.add_node('g1_grid')
        self.add_node('screen_grid')
        self.add_node('plate')
        self.add_node('cathode')
        self.add_node('V_in_AC')

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

        # 7. Resistenza della griglia schermo
        self.add_component(Resistor(name="R_Screen", node1='V_Bplus_node', node2='screen_grid', resistance=RScreen_val))

        # 8. Condensatore di bypass della griglia schermo
        self.add_component(Capacitor(name="C_Screen", node1='screen_grid', node2='gnd', capacitance=CScreen_val))

        # 9. Resistenza anodica
        self.add_component(Resistor(name="R_Plate", node1='V_Bplus_node', node2='plate', resistance=RPlate_val))

        # 10. Pentodo
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

    pentode_exciter = PentodeHarmonicExciter()
    
    print(f"Circuito: {pentode_exciter.name}")
    print("Nodi:", pentode_exciter.get_node_names())
    print("Componenti:")
    for comp in pentode_exciter.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(pentode_exciter)

    print("\nSimulazione DC del Pentode Harmonic Exciter...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in pentode_exciter.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Verifica il punto di riposo del Pentodo
        v1 = next((c for c in pentode_exciter.nonlinear_components if isinstance(c, Pentode) and c.name == "V1"), None)
        if v1:
            V_plate_q = node_voltages[solver._get_node_index(v1.nodes['plate'])]
            V_g1_grid_q = node_voltages[solver._get_node_index(v1.nodes['g1_grid'])]
            V_screen_grid_q = node_voltages[solver._get_node_index(v1.nodes['screen_grid'])]
            V_cathode_q = node_voltages[solver._get_node_index(v1.nodes['cathode'])]
            
            Vgk_q = V_g1_grid_q - V_cathode_q
            Vpk_q = V_plate_q - V_cathode_q
            Vsg_q = V_screen_grid_q - V_cathode_q
            
            Ip_q = v1.calculate_plate_current(Vgk_q, Vpk_q, Vsg_q)

            print(f"\n--- Punto di Riposo del Pentodo (V1) ---")
            print(f"  V_GK (V1) = {Vgk_q:.3f} V")
            print(f"  V_PK (V1) = {Vpk_q:.3f} V")
            print(f"  V_SG (V1) = {Vsg_q:.3f} V")
            print(f"  I_P (V1) = {Ip_q*1e3:.3f} mA")
            
            if Ip_q > 1e-3 and Ip_q < 30e-3: # Corrente anodica per operare in classe A
                print(f"Il Pentodo (V1) è correttamente polarizzato per l'eccitazione armonica.")
            else:
                print(f"ATTENZIONE: Il Pentodo (V1) potrebbe essere fuori dalla polarizzazione desiderata per la distorsione.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
