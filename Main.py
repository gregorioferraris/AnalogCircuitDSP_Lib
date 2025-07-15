import matplotlib.pyplot as plt
import numpy as np

from circuit_solver.circuit import Circuit
from circuit_solver.mna_solver import MnaSolver
from circuit_solver.subcircuit import Subcircuit # Nuovo import

# Importa i componenti
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor
from components.voltage_source import VoltageSource
from components.current_source import CurrentSource
from components.diode import Diode
from components.transformer import Transformer
from components.splitter import Splitter # Nuovo import

# --- Esempio di Subcircuit: Attenuatore (come avevi già fatto) ---
class AttenuatorBlock(Subcircuit):
    def __init__(self, name, input_node_name, output_node_name, ground_node_name, attenuation_factor=0.5, total_resistance=10e3):
        super().__init__(name, input_nodes=[input_node_name], output_nodes=[output_node_name])
        
        r2_val = total_resistance * attenuation_factor
        r1_val = total_resistance - r2_val
        
        if r1_val < 1e-6: r1_val = 1e-6
        if r2_val < 1e-6: r2_val = 1e-6

        # Aggiungi i nodi al sottocircuito
        self.add_node(input_node_name)
        self.add_node(output_node_name)
        self.add_node(ground_node_name) # Assicurati che '0' sia mappato a ground_node_name se diverso

        # Aggiungi i componenti al sottocircuito
        self.add_component(Resistor(name=f"{name}_R1", node1=input_node_name, node2=output_node_name, resistance=r1_val))
        self.add_component(Resistor(name=f"{name}_R2", node1=output_node_name, node2=ground_node_name, resistance=r2_val))
        
        self.input_node_name = input_node_name
        self.output_node_name = output_node_name
        self.ground_node_name = ground_node_name
        self.total_resistance = total_resistance

    def set_attenuation_factor(self, factor):
        if not (0 <= factor <= 1):
            raise ValueError("Attenuation factor must be between 0 and 1.")
        
        r2_val = self.total_resistance * factor
        r1_val = self.total_resistance - r2_val
        if r1_val < 1e-6: r1_val = 1e-6
        if r2_val < 1e-6: r2_val = 1e-6

        # Trova e aggiorna i resistori (assumendo l'ordine di aggiunta)
        for comp in self.components:
            if comp.name == f"{self.name}_R1":
                comp.resistance = r1_val
            elif comp.name == f"{self.name}_R2":
                comp.resistance = r2_val
        print(f"Attenuator '{self.name}' set to factor {factor:.2f}. R1={r1_val:.0f} Ohm, R2={r2_val:.0f} Ohm")

# --- Esempio di Subcircuit: Filtro RC Passa-Basso ---
class LowPassFilterBlock(Subcircuit):
    def __init__(self, name, input_node_name, output_node_name, ground_node_name, R_val, C_val):
        super().__init__(name, input_nodes=[input_node_name], output_nodes=[output_node_name])
        
        self.add_node(input_node_name)
        self.add_node(output_node_name)
        self.add_node(ground_node_name)

        self.add_component(Resistor(name=f"{name}_R", node1=input_node_name, node2=output_node_name, resistance=R_val))
        self.add_component(Capacitor(name=f"{name}_C", node1=output_node_name, node2=ground_node_name, capacitance=C_val))
        
        self.input_node_name = input_node_name
        self.output_node_name = output_node_name
        self.ground_node_name = ground_node_name

# --- Sorgente di Tensione Sinusoidale (già definita nel tuo main.py) ---
class SineVoltageSource(VoltageSource):
    def __init__(self, name, node_plus, node_minus, amplitude, frequency, phase=0.0):
        super().__init__(name, node_plus, node_minus, initial_voltage=0.0)
        self.amplitude = amplitude
        self.frequency = frequency
        self.phase = phase # in radianti
    
    def get_voltage(self, time: float) -> float:
        return self.amplitude * np.sin(2 * np.pi * self.frequency * time + self.phase)


def main():
    # --- Esempio: Catena di Segnale con Attenuatore e Filtro RC ---
    print("\n--- Simulating Signal Chain: Attenuator -> Low-Pass Filter ---")

    main_circuit = Circuit()

    # 1. Crea i blocchi (Subcircuit)
    attenuator_block = AttenuatorBlock("Att_1", "atten_in", "atten_out", "0", attenuation_factor=0.5)
    filter_block = LowPassFilterBlock("LPF_1", "filter_in", "filter_out", "0", R_val=10e3, C_val=10e-9) # 10kOhm, 10nF

    # 2. Connetti i blocchi usando connect_blocks (il "cable" logico)
    # Questo aggiunge i componenti interni dei blocchi al main_circuit e connette i loro nodi I/O
    main_circuit.connect_blocks(attenuator_block, filter_block,
                                source_output_names=["atten_out"],
                                dest_input_names=["filter_in"])

    # 3. Aggiungi una sorgente di ingresso al primo blocco
    main_circuit.add_component(SineVoltageSource("V_in_chain", "atten_in", "0", amplitude=10.0, frequency=500.0)) # 10Vpk, 500 Hz

    # 4. Aggiungi un carico all'uscita finale (opzionale, ma buona pratica)
    main_circuit.add_component(Resistor("R_load", "filter_out", "0", resistance=10e3)) # 10 kOhm

    # Inizializza e risolvi
    solver = MnaSolver(main_circuit)

    sim_time = 0.01 # 10 ms (5 cicli a 500 Hz)
    time_step = 1e-6 # 1 us
    times, solution_history = solver.simulate_transient(0, sim_time, time_step)

    if len(solution_history) > 0:
        # Recupera gli ID dei nodi per il plotting
        input_chain_id = main_circuit.node_map['atten_in']
        output_atten_id = main_circuit.node_map['atten_out'] # Questo è lo stesso di 'filter_in'
        output_filter_id = main_circuit.node_map['filter_out']

        # Estrai le tensioni
        v_input_chain = solution_history[:, input_chain_id]
        v_output_atten = solution_history[:, output_atten_id]
        v_output_filter = solution_history[:, output_filter_id]

        # Plot dei risultati
        plt.figure(figsize=(12, 8))
        plt.plot(times, v_input_chain, label='V_input_chain (V)')
        plt.plot(times, v_output_atten, label='V_output_atten (V) / V_input_filter (V)')
        plt.plot(times, v_output_filter, label='V_output_filter (V)')
        plt.title('Simulazione Catena di Segnale: Attenuatore -> Filtro Passa-Basso')
        plt.xlabel('Tempo (s)')
        plt.ylabel('Tensione (V)')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
    else:
        print("No solution history to plot for signal chain.")

    # Puoi aggiungere qui gli altri esempi (raddrizzatore, RLC) se vuoi mantenerli
    # ... (codice per Half-Wave Rectifier e RLC Circuit dagli esempi precedenti) ...


if __name__ == "__main__":
    main()

