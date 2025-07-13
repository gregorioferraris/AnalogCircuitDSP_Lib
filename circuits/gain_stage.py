# circuits/gain_stage.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.opamp import OpAmp # Richiede la tua classe OpAmp
from circuit_solver.circuit import Circuit
from components.voltage_source import VoltageSource # Per i test

class GainStage(Circuit):
    def __init__(self, name="Gain_Stage", input_node='in', output_node='out', gain=2.0):
        super().__init__(name)

        if gain < 1.0:
            print("Warning: Gain < 1 for a non-inverting amplifier will be ignored, setting to 1.0 (buffer).")
            gain = 1.0 # La configurazione non invertente non può attenuare al di sotto di 1

        self.input_node = input_node
        self.output_node = output_node
        self.gain = gain

        # Nodi
        self.add_node(input_node)
        self.add_node(output_node)
        self.add_node('opamp_in_non_inv') # Ingresso non invertente Op-Amp
        self.add_node('opamp_in_inv')     # Ingresso invertente Op-Amp
        self.add_node('opamp_out_internal') # Uscita interna Op-Amp

        # Valori resistenze per il guadagno (esempio)
        # G = 1 + Rf/Rg
        # Se G = 2, allora Rf/Rg = 1. Possiamo scegliere Rf = 10k, Rg = 10k
        # Se G = 5, allora Rf/Rg = 4. Possiamo scegliere Rf = 40k, Rg = 10k
        
        R_g_val = 10e3 # Resistenza di riferimento per Rg
        if gain == 1.0:
            R_f_val = 1e-6 # Valore molto piccolo, quasi un corto per Rf=0
        else:
            R_f_val = R_g_val * (gain - 1.0)
        
        # Condensatori di accoppiamento
        C_in_val = 1.0e-6
        C_out_val = 1.0e-6

        # --- Componenti ---

        # Condensatore di accoppiamento in ingresso
        self.add_component(Capacitor(name="C_in", node1=input_node, node2='opamp_in_non_inv', capacitance=C_in_val))
        
        # Resistenza di bias per l'ingresso non invertente (a massa, per DC stability)
        self.add_component(Resistor(name="R_bias_non_inv", node1='opamp_in_non_inv', node2='gnd', resistance=1.0e6))

        # Op-Amp
        self.add_component(OpAmp(name="U1_gain",
                                  nodes={'in_inv': 'opamp_in_inv',
                                         'in_non_inv': 'opamp_in_non_inv',
                                         'out': 'opamp_out_internal'}))

        # Rete di feedback per il guadagno
        self.add_component(Resistor(name="R_f", node1='opamp_out_internal', node2='opamp_in_inv', resistance=R_f_val))
        self.add_component(Resistor(name="R_g", node1='opamp_in_inv', node2='gnd', resistance=R_g_val))

        # Condensatore di accoppiamento in uscita
        self.add_component(Capacitor(name="C_out", node1='opamp_out_internal', node2=output_node, capacitance=C_out_val))
        self.add_component(Resistor(name="R_load_out", node1=output_node, node2='gnd', resistance=10e3))


# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver
    import matplotlib.pyplot as plt

    # Test del GainStage in un circuito temporale
    class TestGainCircuit(Circuit):
        def __init__(self, gain_val):
            super().__init__("Test_Gain_Circuit")
            self.add_node('sig_input')
            self.add_node('sig_output')
            
            # Sorgente di segnale AC
            self.input_voltage_source = VoltageSource(name="V_in", nodes={'pos': 'sig_input', 'neg': 'gnd'}, voltage=0.0)
            self.input_voltage_source.get_voltage = lambda t: 0.5 * np.sin(2 * np.pi * 1000 * t) # 1kHz sine wave
            self.add_component(self.input_voltage_source)

            # Stadio di Guadagno
            self.gain_stage = GainStage(name="Amplifier_1", input_node='sig_input', output_node='sig_output', gain=gain_val)
            for comp in self.gain_stage.components:
                self.add_component(comp)
            
            # Assicurati che le alimentazioni OpAmp siano definite nel tuo circuito se il modello non è ideale
            # Per l'OpAmp ideale, non necessitano di Vs aggiuntive, ma i nodi 'VCC_pos', 'VCC_neg' possono essere usati.
            # Se la tua classe OpAmp richiede nodi Vcc, li dovresti aggiungere qui:
            # self.add_component(VoltageSource(name="VCC_Pos_OpAmp", nodes={'pos': 'VCC_pos', 'neg': 'gnd'}, voltage=15.0))
            # self.add_component(VoltageSource(name="VCC_Neg_OpAmp", nodes={'pos': 'VCC_neg', 'neg': 'gnd'}, voltage=-15.0))
            # E modificare la classe OpAmp per usare questi nodi.


    gain_to_test = 2.0 # Prova un guadagno di 2x

    test_circuit = TestGainCircuit(gain_to_test)
    solver = MnaSolver(test_circuit)
    
    print(f"\nSimulazione transitoria del Gain Stage (Gain={gain_to_test}x)...")
    t_start = 0
    t_end = 0.005 # 5 ms, 5 cicli a 1kHz
    dt = 1.0e-6 # 1 microsecondo, 1000 campioni per ciclo

    results, time_points = solver.solve_transient(t_start, t_end, dt, output_nodes=['sig_input', 'sig_output'])

    if results is not None:
        plt.figure(figsize=(12, 6))
        plt.plot(results['time'], results['sig_input'], label='Input Signal (V_sig_input)')
        plt.plot(results['time'], results['sig_output'], label='Output Signal (V_sig_output)')
        plt.title(f'Gain Stage Transient Simulation (Gain={gain_to_test}x)')
        plt.xlabel('Time (s)')
        plt.ylabel('Voltage (V)')
        plt.grid(True)
        plt.legend()
        plt.show()

        # Verifica il guadagno medio dopo un po' di tempo
        # Ignora i primi cicli per la stabilizzazione
        start_idx = int(0.002 / dt)
        avg_gain = np.mean(np.abs(np.array(results['sig_output'][start_idx:]) / np.array(results['sig_input'][start_idx:])))
        print(f"\nGuadagno medio calcolato: {avg_gain:.2f}x (atteso: {gain_to_test:.2f}x)")

    else:
        print("La simulazione transitoria del Gain Stage non è riuscita.")
