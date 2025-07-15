import matplotlib.pyplot as plt
import numpy as np

from circuit_solver.circuit import Circuit
from circuit_solver.mna_solver import MnaSolver
from circuit_solver.subcircuit import Subcircuit

# Importa i componenti di base
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor
from components.voltage_source import VoltageSource

# Importa i blocchi di stadio analogico (Subcircuit)
from components.analog_input_stage import AnalogInputStage
from components.analog_output_stage import AnalogOutputStage

# Importa il modulo di conversione digitale (funzione)
from utils.digital_converter import simulate_vintage_adc_dac

# Importa i nuovi componenti funzionali e il PluginProcessor
from components.external_effect_wrapper import ExternalEffectWrapper
from components.splitter import Splitter # Splitter è ora un componente MNA!
from components.summation import Summation # Rimane funzionale
from utils.plugin_processor import PluginProcessor, ArrayVoltageSource

# Classe helper per sorgente sinusoidale (già definita)
class SineVoltageSource(VoltageSource):
    def __init__(self, name, node_plus, node_minus, amplitude, frequency, phase=0.0):
        super().__init__(name, node_plus, node_minus, initial_voltage=0.0)
        self.amplitude = amplitude
        self.frequency = frequency
        self.phase = phase # in radianti
    
    def get_voltage(self, time: float) -> float:
        return self.amplitude * np.sin(2 * np.pi * self.frequency * time + self.phase)


def main():
    print("\n--- Simulating Circuit with Complex Plugin Workflow ---")

    sim_time = 0.02 # 20 ms
    time_step = 20e-6 # 20 us (simula una "risoluzione" temporale per l'MNA)
    sample_rate_for_processing = int(1.0 / time_step) # Frequenza di campionamento per gli effetti esterni

    # --- DEFINIZIONE DEI BLOCCHI DI PROCESSING E ROUTING ---
    
    # 1. Stadio di Ingresso Analogico (Subcircuit MNA)
    input_analog_stage = AnalogInputStage(
        "PreAmp_Analog_In", "input_signal_raw", "analog_stage_out", "0",
        input_resistance=10e3,
        gain_factor=1.5, # Un po' di guadagno per spingere
        filter_R=100, # Filtro leggero
        filter_C=100e-9,
        diode_saturation_current=1e-12
    )

    # 2. Splitter (ORA È UN COMPONENTE MNA!)
    # Prende l'output dello stadio analogico di ingresso e lo splitta
    # Le sue uscite saranno nodi reali nel circuito MNA.
    main_splitter = Splitter("MainSplitter", "analog_stage_out", "splitter_out_a_node", "splitter_out_b_node")

    # 3. Processing A (Interno alla libreria - Distorsione) - Rimane funzionale
    # Questo effetto prenderà il segnale dal nodo 'splitter_out_a_node'
    distortion_effect_a = ExternalEffectWrapper(
        "Distortion_A", "dist_in_a", "dist_out_a", # Nodi concettuali per l'effetto
        effect_type='distortion_simple',
        param1=8.0, # Fattore di drive più alto
        sample_rate=sample_rate_for_processing
    )

    # 4. Processing B (Interno alla libreria - Ritardo) - Rimane funzionale
    # Questo effetto prenderà il segnale dal nodo 'splitter_out_b_node'
    delay_effect_b = ExternalEffectWrapper(
        "Delay_B", "delay_in_b", "delay_out_b", # Nodi concettuali per l'effetto
        effect_type='delay',
        param1=0.005, # 5 ms di ritardo
        sample_rate=sample_rate_for_processing
    )
    delay_effect_b.reset_state() 

    # 5. Summation (Funzionale) - Rimane funzionale
    # Somma l'output del processing A e l'output del processing B
    main_summation = Summation("MainSummation", ["sum_in_1", "sum_in_2"], "sum_out")

    # 6. Stadio di Uscita Analogico (Subcircuit MNA)
    output_analog_stage = AnalogOutputStage(
        "PostAmp_Analog_Out", "final_analog_in", "output_final", "0",
        filter_R=500,
        filter_C=220e-9,
        output_resistance=50.0,
        clipper_threshold=0.7, # Clipping più leggero
        diode_saturation_current=1e-9
    )

    # --- CONFIGURAZIONE E ESECUZIONE DEL FLUSSO DI LAVORO ---

    # 1. Ingresso Raw -> Stadio Analogico di Ingresso (MNA) -> Splitter (MNA)
    print("\n-- Fase 1: Ingresso Raw -> Stadio Analogico di Ingresso -> Splitter (Tutto MNA) --")
    circuit_segment_1 = Circuit()
    circuit_segment_1.add_block(input_analog_stage)
    circuit_segment_1.add_component(SineVoltageSource("Audio_In", "input_signal_raw", "0", amplitude=0.8, frequency=200.0)) # Segnale più moderato
    
    # Aggiungi lo Splitter direttamente al circuito MNA
    circuit_segment_1.add_component(main_splitter)

    solver_1 = MnaSolver(circuit_segment_1)
    times_segment_1, solution_history_1 = solver_1.simulate_transient(0, sim_time, time_step)
    
    # Ora i segnali splittati vengono letti direttamente dai nodi di uscita dello Splitter nel circuito MNA
    signal_splitter_out_a = solution_history_1[:, circuit_segment_1.node_map['splitter_out_a_node']]
    signal_splitter_out_b = solution_history_1[:, circuit_segment_1.node_map['splitter_out_b_node']]
    
    # 2. Splitter Output A -> S Processing A (PluginProcessor -> ExternalEffectWrapper)
    # Il segnale 'signal_splitter_out_a' è il "cable" che entra nel primo processore funzionale
    print("\n-- Fase 2: Splitter Output A -> S Processing A (PluginProcessor) --")
    processor_a = PluginProcessor(distortion_effect_a, sample_rate_for_processing)
    processed_signal_a_vs = processor_a.process_segment(signal_splitter_out_a, times_segment_1)
    
    # 3. Splitter Output B -> S Processing B (PluginProcessor -> ExternalEffectWrapper)
    # Il segnale 'signal_splitter_out_b' è il "cable" che entra nel secondo processore funzionale
    print("\n-- Fase 3: Splitter Output B -> S Processing B (PluginProcessor) --")
    processor_b = PluginProcessor(delay_effect_b, sample_rate_for_processing)
    processed_signal_b_vs = processor_b.process_segment(signal_splitter_out_b, times_segment_1)

    # 4. Cable -> Summation (Funzionale)
    # I segnali elaborati (voltages dagli ArrayVoltageSource) sono i nostri "cable"
    print("\n-- Fase 4: Cable -> Summation (Funzionale) --")
    summation_inputs = {
        "sum_in_1": processed_signal_a_vs.voltages,
        "sum_in_2": processed_signal_b_vs.voltages
    }
    final_summed_signal = main_summation.process_block(summation_inputs)

    # 5. Cable -> Stadio Analogico di Uscita (MNA)
    print("\n-- Fase 5: Cable -> Stadio Analogico di Uscita (MNA) --")
    circuit_segment_2 = Circuit()
    circuit_segment_2.add_block(output_analog_stage)
    
    # Inietta il segnale sommato come input per lo stadio finale MNA
    circuit_segment_2.add_component(ArrayVoltageSource(
        "Summed_Signal_VS",
        output_analog_stage.input_nodes[0], # Connetti all'input dello stadio di uscita
        "0",
        times=times_segment_1, # Usa gli stessi tempi della simulazione originale
        voltages=final_summed_signal
    ))
    circuit_segment_2.add_component(Resistor("Output_Load", "output_final", "0", resistance=10e3))

    solver_2 = MnaSolver(circuit_segment_2)
    times_segment_2, solution_history_2 = solver_2.simulate_transient(0, sim_time, time_step)
    
    # --- ESTRAZIONE E PLOTTING DEI RISULTATI ---
    print("\n-- Fase 6: Plotting dei Risultati --")
    if len(solution_history_2) > 0:
        # Segnali per il plotting
        v_input_raw = solution_history_1[:, circuit_segment_1.node_map['input_signal_raw']]
        v_pre_splitter = solution_history_1[:, circuit_segment_1.node_map['analog_stage_out']] # Output dello stadio analogico in ingresso
        v_splitter_out_a = solution_history_1[:, circuit_segment_1.node_map['splitter_out_a_node']] # Output dello splitter A
        v_splitter_out_b = solution_history_1[:, circuit_segment_1.node_map['splitter_out_b_node']] # Output dello splitter B

        v_post_distortion_a = processed_signal_a_vs.voltages
        v_post_delay_b = processed_signal_b_vs.voltages
        v_summed = final_summed_signal
        v_final_output = solution_history_2[:, circuit_segment_2.node_map['output_final']]

        plt.figure(figsize=(14, 12))
        
        plt.subplot(4, 1, 1)
        plt.plot(times_segment_1, v_input_raw, label='1. Ingresso Raw (V)', alpha=0.7)
        plt.plot(times_segment_1, v_pre_splitter, label='2. Uscita Pre-Splitter (V)', linestyle='--', alpha=0.8)
        plt.title('Ingresso e Uscita Stadio Analogico Iniziale (Pre-Splitter)')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensione (V)')
        plt.grid(True)
        plt.legend()

        plt.subplot(4, 1, 2)
        plt.plot(times_segment_1, v_pre_splitter, label='2. Ingresso Splitter (V)', linestyle='--', alpha=0.8)
        plt.plot(times_segment_1, v_splitter_out_a, label='2a. Uscita Splitter A (V)', linestyle=':', linewidth=2, color='orange')
        plt.plot(times_segment_1, v_splitter_out_b, label='2b. Uscita Splitter B (V)', linestyle='-.', linewidth=2, color='green')
        plt.title('Comportamento dello Splitter MNA')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensione (V)')
        plt.grid(True)
        plt.legend()

        plt.subplot(4, 1, 3)
        plt.plot(times_segment_1, v_splitter_out_a, label='3a. Ingresso Distorsione (V)', linestyle=':', linewidth=2, color='orange')
        plt.plot(times_segment_1, v_post_distortion_a, label='3a. Uscita Distorsione (V)', linestyle='--', alpha=0.8, color='red')
        plt.plot(times_segment_1, v_splitter_out_b, label='3b. Ingresso Ritardo (V)', linestyle='-.', linewidth=2, color='green')
        plt.plot(times_segment_1, v_post_delay_b, label='3b. Uscita Ritardo (V)', linestyle='--', alpha=0.8, color='blue')
        plt.title('Elaborazione dei Percorsi Paralleli')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensione (V)')
        plt.grid(True)
        plt.legend()

        plt.subplot(4, 1, 4)
        plt.plot(times_segment_1, v_post_distortion_a, label='4a. Ingresso Sommatoria (Distorsione) (V)', linestyle=':', linewidth=2, color='orange')
        plt.plot(times_segment_1, v_post_delay_b, label='4b. Ingresso Sommatoria (Ritardo) (V)', linestyle='-.', linewidth=2, color='green')
        plt.plot(times_segment_1, v_summed, label='5. Segnale Sommato (V)', linewidth=2, color='blue')
        plt.plot(times_segment_2, v_final_output, label='6. Uscita Finale (Post-Analog Stage) (V)', linewidth=2, color='red')
        plt.title('Sommatoria e Uscita Finale')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensione (V)')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()
    else:
        print("Nessuna storia della soluzione da plottare.")

if __name__ == "__main__":
    main()

