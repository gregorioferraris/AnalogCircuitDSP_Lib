# circuits/triode_bass_enhancer.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.triode import Triode # Assicurati che la tua classe Triode sia importabile
from circuit_solver.circuit import Circuit

class TriodeBassEnhancer(Circuit):
    def __init__(self, name="Triode_Bass_Enhancer"):
        super().__init__(name)

        # --- Parametri del circuito ---
        # 1. Filtro Passa-Basso (LPF) RC per isolare le basse frequenze
        # Frequenza di taglio circa 1/(2*pi*R_lpf*C_lpf)
        # R_lpf1 = 10k, C_lpf1 = 0.47uF -> Fc ~ 33 Hz (per i bassi molto profondi)
        # R_lpf2 = 10k, C_lpf2 = 0.1uF  -> Fc ~ 159 Hz (per bassi più "punchy")
        R_lpf1_val = 10e3
        C_lpf1_val = 0.47e-6
        R_lpf2_val = 10e3
        C_lpf2_val = 0.1e-6 # Questo rende il filtro un po' più aggressivo sui medi-bassi

        # 2. Stadio a Triodo per Saturazione
        V_Bplus_val = 150.0  # Alimentazione anodica (più bassa per saturare più facilmente)
        RG_val = 1.0e6       # Resistenza di griglia
        RK_val = 1.5e3       # Resistenza di catodo (bias)
        RPlate_val = 68e3    # Resistenza anodica (carico)

        # Parametri del Triodo (es. 12AU7 o 12AT7 - triodi a guadagno medio/basso, più adatti alla saturazione morbida)
        triode_params = {
            "mu": 17.0,        # Fattore di amplificazione (es. 12AU7)
            "rp": 7.7e3,       # Resistenza interna della placca (Ohm, 12AU7)
            "gm": 0.0022      # Transconduttanza (Siemens, 12AU7)
        }
        
        # Condensatori di accoppiamento
        C_in_val = 0.1e-6    # Ingresso LPF
        C_triode_in_val = 0.1e-6 # Ingresso Triodo
        C_triode_out_val = 0.1e-6 # Uscita Triodo
        CK_val = 47e-6       # Bypass catodo (più grande per massimo guadagno anche sui bassi)

        # 3. Sezione Mixer (molto semplificata con resistori)
        # Per un mixer più robusto servirebbero Op-Amp, ma per simulazione possiamo usare nodi separati e resistori.
        R_mix_dry_val = 10e3   # Resistenza per il percorso Dry
        R_mix_wet_val = 10e3   # Resistenza per il percorso Wet (saturato)
        R_out_load_val = 10e3  # Resistenza di carico in uscita

        # --- Nodi del circuito ---
        self.add_node('in')             # Ingresso segnale AC
        self.add_node('lpf_out')        # Uscita filtro passa-basso (segnale per il triodo)
        self.add_node('V_Bplus_node')   # Alimentazione triodo
        self.add_node('triode_grid')    # Griglia del triodo
        self.add_node('triode_plate')   # Anodo del triodo
        self.add_node('triode_cathode') # Catodo del triodo
        self.add_node('triode_out')     # Uscita AC del triodo (segnale wet)
        self.add_node('mixer_node')     # Nodo di somma del mixer
        self.add_node('out')            # Uscita finale del basso enhancer
        self.add_node('V_in_AC')        # Nodo per la sorgente di segnale AC

        # --- Componenti ---

        # Sorgente di segnale d'ingresso AC
        self.add_component(VoltageSource(name="V_signal", nodes={'pos': 'V_in_AC', 'neg': 'gnd'}, voltage=0.0))
        self.add_component(Capacitor(name="C_in_main", node1='V_in_AC', node2='in', capacitance=C_in_val)) # Accoppiamento ingresso principale

        # --- 1. Sezione Filtro Passa-Basso (LPF) ---
        # Filtro RC passa-basso per isolare le basse frequenze
        self.add_component(Resistor(name="R_lpf1", node1='in', node2='lpf_node_mid', resistance=R_lpf1_val))
        self.add_component(Capacitor(name="C_lpf1", node1='lpf_node_mid', node2='gnd', capacitance=C_lpf1_val))
        self.add_component(Resistor(name="R_lpf2", node1='lpf_node_mid', node2='lpf_out', resistance=R_lpf2_val))
        self.add_component(Capacitor(name="C_lpf2", node1='lpf_out', node2='gnd', capacitance=C_lpf2_val))


        # --- 2. Sezione Stadio a Triodo per Saturazione ---
        self.add_component(VoltageSource(name="V_Bplus", nodes={'pos': 'V_Bplus_node', 'neg': 'gnd'}, voltage=V_Bplus_val))

        # Accoppiamento ingresso triodo dal LPF
        self.add_component(Capacitor(name="C_triode_in", node1='lpf_out', node2='triode_grid', capacitance=C_triode_in_val))
        self.add_component(Resistor(name="RG_triode", node1='triode_grid', node2='gnd', resistance=RG_val)) # Resistenza di griglia per bias

        # Resistenza di catodo per auto-polarizzazione
        self.add_component(Resistor(name="R_K_triode", node1='triode_cathode', node2='gnd', resistance=RK_val))
        # Condensatore di bypass del catodo
        self.add_component(Capacitor(name="C_K_triode", node1='triode_cathode', node2='gnd', capacitance=CK_val))

        # Resistenza anodica (Plate Resistor)
        self.add_component(Resistor(name="R_Plate_triode", node1='V_Bplus_node', node2='triode_plate', resistance=RPlate_val))

        # Il Triodo stesso
        self.add_component(Triode(name="V1_saturator",
                                  nodes={'plate': 'triode_plate', 'grid': 'triode_grid', 'cathode': 'triode_cathode'},
                                  **triode_params))

        # Condensatore di accoppiamento in uscita dal triodo
        self.add_component(Capacitor(name="C_triode_out", node1='triode_plate', node2='triode_out', capacitance=C_triode_out_val))
        
        # Resistenza di carico per lo stadio di saturazione del triodo
        self.add_component(Resistor(name="R_load_triode", node1='triode_out', node2='gnd', resistance=100e3))

        # --- 3. Sezione Mixer Semplificata ---
        # Percorso DRY (segnale originale)
        self.add_component(Resistor(name="R_mix_dry", node1='in', node2='mixer_node', resistance=R_mix_dry_val))
        
        # Percorso WET (segnale saturato dal triodo)
        self.add_component(Resistor(name="R_mix_wet", node1='triode_out', node2='mixer_node', resistance=R_mix_wet_val))
        
        # Resistenza di carico sul nodo mixer
        self.add_component(Resistor(name="R_mixer_load", node1='mixer_node', node2='gnd', resistance=R_out_load_val))
        
        # Uscita finale del circuito
        self.add_component(Resistor(name="R_final_out", node1='mixer_node', node2='out', resistance=10)) # Piccola resistenza per l'uscita
        self.add_component(Resistor(name="R_final_load_global", node1='out', node2='gnd', resistance=10e3)) # Carico generale


# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    bass_enhancer_valve = TriodeBassEnhancer()
    
    print(f"Circuito: {bass_enhancer_valve.name}")
    print("Nodi:", bass_enhancer_valve.get_node_names())
    print("Componenti:")
    for comp in bass_enhancer_valve.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(bass_enhancer_valve)

    print("\nSimulazione DC del Triode Bass Enhancer...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in bass_enhancer_valve.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Verifica il punto di riposo del Triodo
        v1 = next((c for c in bass_enhancer_valve.nonlinear_components if isinstance(c, Triode) and c.name == "V1_saturator"), None)
        if v1:
            V_plate_q = node_voltages[solver._get_node_index(v1.nodes['triode_plate'])]
            V_grid_q = node_voltages[solver._get_node_index(v1.nodes['triode_grid'])]
            V_cathode_q = node_voltages[solver._get_node_index(v1.nodes['triode_cathode'])]
            
            Vgk_q = V_grid_q - V_cathode_q
            Vpk_q = V_plate_q - V_cathode_q
            Ip_q = v1.calculate_plate_current(Vgk_q, Vpk_q)

            print(f"\n--- Punto di Riposo del Triodo (V1_saturator) ---")
            print(f"  V_GK (V1) = {Vgk_q:.3f} V")
            print(f"  V_PK (V1) = {Vpk_q:.3f} V")
            print(f"  I_P (V1) = {Ip_q*1e3:.3f} mA")
            
            # Condizioni per il funzionamento in regione attiva e per la distorsione
            # Vgk_q dovrebbe essere leggermente negativo, Vpk_q positivo e significativo.
            if Ip_q > 0.5e-3 and Ip_q < 10e-3 and Vgk_q < 0 and Vpk_q > 50:
                print(f"Il Triodo (V1_saturator) è correttamente polarizzato per saturazione di classe A.")
            else:
                print(f"ATTENZIONE: Il Triodo (V1_saturator) potrebbe essere fuori dalla polarizzazione desiderata per la distorsione.")
        
        # Verifica che il nodo del mixer sia vicino a 0V in DC (prima di segnale AC)
        v_mixer_node_dc = node_voltages[solver._get_node_index('mixer_node')]
        print(f"\nV(mixer_node) DC: {v_mixer_node_dc:.3f} V")
        if abs(v_mixer_node_dc) < 1e-6:
            print("Il nodo del mixer è correttamente a 0V in DC.")
        else:
            print("ATTENZIONE: Il nodo del mixer non è a 0V in DC, potrebbe esserci un problema di bias.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
