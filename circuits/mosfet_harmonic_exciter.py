# circuits/mosfet_harmonic_exciter.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.mosfet import MOSFET # Assicurati che la tua classe MOSFET sia importabile
from circuit_solver.circuit import Circuit

class MosfetHarmonicExciter(Circuit):
    def __init__(self, name="MOSFET_Harmonic_Exciter"):
        super().__init__(name)

        # --- Parametri del circuito ---
        VCC_val = 12.0     # 12V - Tensione di alimentazione
        RD_val = 10e3      # 10 kOhm - Resistenza di Drain
        RG_val = 1.0e6     # 1 MOhm - Resistenza di Griglia (alta impedenza d'ingresso)
        RS_val = 1.0e3     # 1 kOhm - Resistenza di Source (per auto-polarizzazione e stabilizzazione)

        C_in_val = 1.0e-6  # 1 uF - Condensatore di accoppiamento in ingresso
        C_out_val = 1.0e-6 # 1 uF - Condensatore di accoppiamento in uscita
        CS_val = 100e-6    # 100 uF - Condensatore di bypass del Source (per massimo guadagno AC)

        # Parametri del MOSFET (es. BS170 o simili)
        mosfet_params = {
            "k": 0.003,     # 3 mA/V^2 - Parametro di transconduttanza
            "Vth": 1.5,     # 1.5 V - Tensione di soglia
            "lambda_val": 0.02 # 0.02 V^-1 - Parametro di modulazione della lunghezza del canale
        }

        # --- Nodi del circuito ---
        self.add_node('in')          # Ingresso del segnale AC
        self.add_node('out')         # Uscita del segnale processato
        self.add_node('VCC_node')    # Nodo di alimentazione positiva
        self.add_node('drain')       # Nodo Drain del MOSFET
        self.add_node('gate')        # Nodo Gate del MOSFET
        self.add_node('source')      # Nodo Source del MOSFET
        self.add_node('V_in_AC')     # Nodo per la sorgente di segnale AC

        # --- Componenti ---

        # 1. Alimentazione DC
        self.add_component(VoltageSource(name="VCC_Power", nodes={'pos': 'VCC_node', 'neg': 'gnd'}, voltage=VCC_val))

        # 2. Sorgente di segnale d'ingresso (per analisi AC/Transitoria)
        self.add_component(VoltageSource(name="V_signal", nodes={'pos': 'V_in_AC', 'neg': 'gnd'}, voltage=0.0))

        # 3. Condensatore di accoppiamento in ingresso
        self.add_component(Capacitor(name="C_in", node1='V_in_AC', node2='in', capacitance=C_in_val))

        # 4. Resistenza di Griglia (polarizzazione Gate a massa per DC)
        self.add_component(Resistor(name="RG", node1='in', node2='gate', resistance=RG_val)) # O direttamente 'in' a 'gate' per isolare DC

        # 5. Resistenza di Drain (collega Drain a VCC)
        self.add_component(Resistor(name="RD", node1='VCC_node', node2='drain', resistance=RD_val))

        # 6. Resistenza di Source (collega Source a gnd, per auto-polarizzazione)
        self.add_component(Resistor(name="RS", node1='source', node2='gnd', resistance=RS_val))

        # 7. Condensatore di bypass sul Source (per alto guadagno AC)
        self.add_component(Capacitor(name="CS", node1='source', node2='gnd', capacitance=CS_val))
        
        # 8. MOSFET (N-Channel)
        self.add_component(MOSFET(name="M1",
                                  nodes={'drain': 'drain', 'gate': 'gate', 'source': 'source'},
                                  **mosfet_params))

        # 9. Condensatore di accoppiamento in uscita (blocca DC, passa AC)
        self.add_component(Capacitor(name="C_out", node1='drain', node2='out', capacitance=C_out_val))
        
        # Aggiungi una resistenza di carico in uscita, se non c'è un carico esterno previsto
        self.add_component(Resistor(name="R_load_out", node1='out', node2='gnd', resistance=100e3)) # 100k Ohm di carico

# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    exciter_circuit = MosfetHarmonicExciter()
    
    print(f"Circuito: {exciter_circuit.name}")
    print("Nodi:", exciter_circuit.get_node_names())
    print("Componenti:")
    for comp in exciter_circuit.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(exciter_circuit)

    print("\nSimulazione DC del MOSFET Harmonic Exciter...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in exciter_circuit.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # Verifica il punto di riposo del MOSFET
        m1 = next((c for c in exciter_circuit.nonlinear_components if isinstance(c, MOSFET) and c.name == "M1"), None)
        if m1:
            V_drain_q = node_voltages[solver._get_node_index(m1.nodes['drain'])]
            V_gate_q = node_voltages[solver._get_node_index(m1.nodes['gate'])]
            V_source_q = node_voltages[solver._get_node_index(m1.nodes['source'])]
            
            Vgs_q = V_gate_q - V_source_q
            Vds_q = V_drain_q - V_source_q
            Id_q = m1.calculate_drain_current(Vgs_q, Vds_q)

            print(f"\n--- Punto di Riposo del MOSFET (M1) ---")
            print(f"  V_GS (M1) = {Vgs_q:.3f} V")
            print(f"  V_DS (M1) = {Vds_q:.3f} V")
            print(f"  I_D (M1) = {Id_q*1e3:.3f} mA")
            
            # Condizione per regione di saturazione (buona per amplificazione lineare, ma qui vogliamo non-linearità)
            if Vgs_q > m1.Vth and Vds_q >= (Vgs_q - m1.Vth):
                print(f"  M1 è in regione di Saturazione (per l'amplificazione).")
            elif Vgs_q <= m1.Vth:
                print(f"  M1 è in Cut-off (spento).")
            else:
                print(f"  M1 è in Triodo (Lineare).")

            if Id_q > 0.05e-3: # Almeno 0.05mA per considerarlo attivo
                print(f"Il MOSFET (M1) è correttamente polarizzato per operare in classe A/AB.")
            else:
                print(f"ATTENZIONE: Il MOSFET (M1) potrebbe essere in cutoff o quasi spento.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
