# circuits/bass_enhancer.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.mosfet import MOSFET # O Diode, Triode, Pentode
from circuit_solver.circuit import Circuit

class BassEnhancer(Circuit):
    def __init__(self, name="Bass_Enhancer"):
        super().__init__(name)

        # --- Parametri del circuito ---
        # Filtro Passa-Basso (LPF) RC per isolare le basse frequenze
        R_lpf = 10e3        # Resistenza LPF
        C_lpf = 0.47e-6     # Condensatore LPF (F_cut = 1 / (2*pi*R*C) ~ 33 Hz)

        # Stadio di Saturazione (Usiamo un MOSFET per la sua non linearità "calda")
        VCC_val = 9.0       # Alimentazione per il MOSFET (più bassa per saturare più facilmente)
        RD_val = 4.7e3      # Resistenza di Drain (valore critico per la polarizzazione)
        RG_val = 1.0e6      # Resistenza di Griglia
        RS_val = 1.0e3      # Resistenza di Source
        CS_val = 100e-6     # Condensatore di bypass del Source

        mosfet_params = {
            "k": 0.003,     # Transconduttanza
            "Vth": 1.5,     # Tensione di soglia
            "lambda_val": 0.02 # Modulazione lunghezza canale
        }
        
        # Resistenza di carico per il MOSFET
        R_load_mosfet = 22e3 # 22 kOhm

        # Condensatori di accoppiamento
        C_in_val = 1.0e-6
        C_out_val = 1.0e-6

        # --- Nodi del circuito ---
        self.add_node('in')             # Ingresso segnale AC
        self.add_node('lpf_out')        # Uscita filtro passa-basso
        self.add_node('VCC_node')       # Alimentazione MOSFET
        self.add_node('drain')          # Drain MOSFET
        self.add_node('gate')           # Gate MOSFET
        self.add_node('source')         # Source MOSFET
        self.add_node('saturator_out')  # Uscita dello stadio di saturazione
        self.add_node('out')            # Uscita finale (dopo un eventuale mixer, non incluso qui)
        self.add_node('V_in_AC')        # Nodo per la sorgente di segnale AC

        # --- Componenti ---

        # Sorgente di segnale d'ingresso AC
        self.add_component(VoltageSource(name="V_signal", nodes={'pos': 'V_in_AC', 'neg': 'gnd'}, voltage=0.0))
        self.add_component(Resistor(name="R_in_sig", node1='V_in_AC', node2='in', resistance=1e3))

        # 1. Filtro Passa-Basso (LPF)
        # Il segnale "in" viene filtrato per isolare solo le basse frequenze.
        self.add_component(Resistor(name="R_lpf", node1='in', node2='lpf_out', resistance=R_lpf))
        self.add_component(Capacitor(name="C_lpf", node1='lpf_out', node2='gnd', capacitance=C_lpf))

        # 2. Alimentazione DC per lo stadio di saturazione
        self.add_component(VoltageSource(name="VCC_Power", nodes={'pos': 'VCC_node', 'neg': 'gnd'}, voltage=VCC_val))

        # 3. Stadio di Saturazione (basato su MOSFET)
        # Condensatore di accoppiamento per il gate del MOSFET
        self.add_component(Capacitor(name="C_in_mosfet", node1='lpf_out', node2='gate', capacitance=C_in_val))
        self.add_component(Resistor(name="RG", node1='gate', node2='gnd', resistance=RG_val)) # Pull-down per bias

        # Resistenza di Drain
        self.add_component(Resistor(name="RD", node1='VCC_node', node2='drain', resistance=RD_val))
        # Resistenza di Source
        self.add_component(Resistor(name="RS", node1='source', node2='gnd', resistance=RS_val))
        # Condensatore di bypass del Source
        self.add_component(Capacitor(name="CS", node1='source', node2='gnd', capacitance=CS_val))
        
        # Il MOSFET stesso
        self.add_component(MOSFET(name="M1_saturator",
                                  nodes={'drain': 'drain', 'gate': 'gate', 'source': 'source'},
                                  **mosfet_params))

        # Condensatore di accoppiamento in uscita dal MOSFET
        self.add_component(Capacitor(name="C_out_mosfet", node1='drain', node2='saturator_out', capacitance=C_out_val))
        
        # Resistenza di carico per lo stadio di saturazione
        self.add_component(Resistor(name="R_load_saturator", node1='saturator_out', node2='gnd', resistance=R_load_mosfet))

        # Output del Basso Enhancer (potrebbe andare a un mixer esterno con il segnale dry)
        # Qui l'uscita è semplicemente il segnale saturato
        self.add_component(Resistor(name="R_final_out", node1='saturator_out', node2='out', resistance=1e3))
        self.add_component(Resistor(name="R_final_load", node1='out', node2='gnd', resistance=10e3))


# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    bass_enhancer = BassEnhancer()
    
    print(f"Circuito: {bass_enhancer.name}")
    print("Nodi:", bass_enhancer.get_node_names())
    print("Componenti:")
    for comp in bass_enhancer.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(bass_enhancer)

    print("\nSimulazione DC del Bass Enhancer...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in bass_enhancer.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Verifica il punto di riposo del MOSFET
        m1 = next((c for c in bass_enhancer.nonlinear_components if isinstance(c, MOSFET) and c.name == "M1_saturator"), None)
        if m1:
            V_drain_q = node_voltages[solver._get_node_index(m1.nodes['drain'])]
            V_gate_q = node_voltages[solver._get_node_index(m1.nodes['gate'])]
            V_source_q = node_voltages[solver._get_node_index(m1.nodes['source'])]
            
            Vgs_q = V_gate_q - V_source_q
            Vds_q = V_drain_q - V_source_q
            Id_q = m1.calculate_drain_current(Vgs_q, Vds_q)

            print(f"\n--- Punto di Riposo del MOSFET (M1_saturator) ---")
            print(f"  V_GS (M1) = {Vgs_q:.3f} V")
            print(f"  V_DS (M1) = {Vds_q:.3f} V")
            print(f"  I_D (M1) = {Id_q*1e3:.3f} mA")
            
            if Id_q > 0.05e-3 and Vgs_q > m1.Vth : # Almeno 0.05mA e sopra la soglia
                print(f"Il MOSFET (M1_saturator) è correttamente polarizzato per saturazione.")
            else:
                print(f"ATTENZIONE: Il MOSFET (M1_saturator) potrebbe essere in cutoff o non saturare.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
