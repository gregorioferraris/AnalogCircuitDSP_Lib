# circuits/transient_detector.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.diode import Diode
from components.voltage_source import VoltageSource # Per il segnale di ingresso
from circuit_solver.circuit import Circuit
# Avresti bisogno di una classe OpAmp per il comparatore, ma non l'abbiamo ancora.
# Qui simuleremo solo la parte analogica di inviluppo e differenziazione.

class TransientDetector(Circuit):
    def __init__(self, name="Transient_Detector"):
        super().__init__(name)

        # --- Parametri del circuito ---
        # Componenti del raddrizzatore a piena onda
        R_in_rect = 10e3 # Resistenza di ingresso al raddrizzatore
        D_params = {"Is": 1.0e-14, "Vt": 0.02585, "N": 1.0} # Parametri diodo generico

        # Componenti del follower di inviluppo (RC LPF)
        R_env = 100e3      # Resistenza filtro inviluppo (carica C_env)
        C_env = 1.0e-6     # Condensatore filtro inviluppo (scarica) - la sua scarica determina l'attack/decay
        
        # Componenti del differenziatore
        C_diff = 0.01e-6   # Condensatore differenziatore
        R_diff = 10e3      # Resistenza differenziatore

        # --- Nodi del circuito ---
        self.add_node('in')             # Ingresso segnale AC
        self.add_node('rect_pos')       # Dopo il primo raddrizzatore
        self.add_node('rect_neg')       # Dopo il secondo raddrizzatore
        self.add_node('rect_out')       # Uscita raddrizzata (inviluppo prima del filtro)
        self.add_node('envelope_out')   # Uscita dell'inviluppo liscio
        self.add_node('transient_out')  # Uscita del rilevatore di transienti (impulso)
        self.add_node('V_in_AC')        # Nodo per la sorgente di segnale AC

        # --- Componenti ---

        # Sorgente di segnale d'ingresso AC
        self.add_component(VoltageSource(name="V_signal", nodes={'pos': 'V_in_AC', 'neg': 'gnd'}, voltage=0.0))
        self.add_component(Resistor(name="R_in_sig", node1='V_in_AC', node2='in', resistance=1e3)) # Resistenza isolante

        # --- Sezione Raddrizzatore a Piena Onda (Bridge Rectifier semplificato) ---
        # Questo è un raddrizzatore a ponte di diodi, per ottenere il valore assoluto
        # Questa è una semplificazione, un raddrizzatore di precisione con Op-Amp sarebbe migliore.
        # D1: Raddrizza il lato positivo del segnale
        self.add_component(Diode(name="D_rect1", nodes={'anode': 'in', 'cathode': 'rect_pos'}, **D_params))
        # D2: Raddrizza il lato negativo del segnale (verso massa)
        self.add_component(Diode(name="D_rect2", nodes={'anode': 'gnd', 'cathode': 'rect_neg'}, **D_params))
        # D3: Connessione per il lato negativo (anodo dal segnale negativo)
        self.add_component(Diode(name="D_rect3", nodes={'anode': 'rect_neg', 'cathode': 'rect_out'}, **D_params))
        # D4: Connessione per il lato positivo (anodo dal segnale positivo)
        self.add_component(Diode(name="D_rect4", nodes={'anode': 'rect_pos', 'cathode': 'rect_out'}, **D_params))
        
        # Pull-down sulla rect_out per riferimento
        self.add_component(Resistor(name="R_rect_pull", node1='rect_out', node2='gnd', resistance=100e3))


        # --- Sezione Segui-Inviluppo (Low-Pass Filter) ---
        # Il condensatore C_env si carica velocemente attraverso la rete del raddrizzatore e si scarica attraverso R_env.
        self.add_component(Resistor(name="R_env", node1='rect_out', node2='envelope_out', resistance=R_env))
        self.add_component(Capacitor(name="C_env", node1='envelope_out', node2='gnd', capacitance=C_env))

        # --- Sezione Differenziatore ---
        # Crea un picco quando l'inviluppo cambia rapidamente
        self.add_component(Capacitor(name="C_diff", node1='envelope_out', node2='transient_out', capacitance=C_diff))
        self.add_component(Resistor(name="R_diff", node1='transient_out', node2='gnd', resistance=R_diff))
        
        # Aggiungiamo una resistenza di carico all'uscita per il tester
        self.add_component(Resistor(name="R_load_transient", node1='transient_out', node2='gnd', resistance=10e3))

# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    transient_detector = TransientDetector()
    
    print(f"Circuito: {transient_detector.name}")
    print("Nodi:", transient_detector.get_node_names())
    print("Componenti:")
    for comp in transient_detector.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(transient_detector)

    print("\nSimulazione DC del Transient Detector...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in transient_detector.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        # In DC, con V_signal a 0V, tutto dovrebbe essere a 0V, e i diodi spenti.
        print(f"\nV(rect_out) DC: {node_voltages[solver._get_node_index('rect_out')]:.3f} V")
        print(f"V(envelope_out) DC: {node_voltages[solver._get_node_index('envelope_out')]:.3f} V")
        print(f"V(transient_out) DC: {node_voltages[solver._get_node_index('transient_out')]:.3f} V")
        
        # Se i valori sono vicini a zero, il DC bias è corretto per il funzionamento AC.
        if all(abs(node_voltages[solver._get_node_index(n)]) < 1e-6 for n in ['rect_out', 'envelope_out', 'transient_out']):
            print("Punto di riposo DC corretto: tutte le uscite sono a zero.")
        else:
            print("ATTENZIONE: Il punto di riposo DC non è a zero, controllare i diodi o le resistenze.")
    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
