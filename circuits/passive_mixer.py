# circuits/passive_mixer.py

from components.resistor import Resistor
from components.capacitor import Capacitor # Per accoppiamento DC se necessario
from components.voltage_source import VoltageSource # Per test
from circuit_solver.circuit import Circuit

class PassiveMixer(Circuit):
    """
    Un sommatore/mixer passivo basato su resistori.
    Semplice attenuazione e somma.
    """
    def __init__(self, name="Passive_Mixer", num_inputs=2, mix_res_val=10e3):
        super().__init__(name)
        
        if num_inputs < 1:
            raise ValueError("Number of inputs must be at least 1.")

        self.num_inputs = num_inputs
        self.mix_res_val = mix_res_val
        self.input_nodes = []

        # Creazione dei nodi di ingresso
        for i in range(num_inputs):
            input_node_name = f'in{i+1}'
            self.add_node(input_node_name)
            self.input_nodes.append(input_node_name)
            
            # Aggiungi una resistenza di miscelazione per ogni ingresso
            self.add_component(Resistor(name=f"R_mix_{i+1}", node1=input_node_name, node2='sum_node', resistance=mix_res_val))

        self.add_node('sum_node') # Nodo comune di somma
        self.add_node('out')      # Uscita AC del mixer

        # Aggiungi una resistenza di carico al nodo di somma per garantire un riferimento.
        # Questo è importante per evitare che il nodo di somma "fluttui" e per stabilire l'impedenza di uscita.
        self.add_component(Resistor(name="R_sum_load", node1='sum_node', node2='gnd', resistance=mix_res_val)) # Stessa resistenza dei mix per semplicità

        # Condensatore di accoppiamento in uscita (per bloccare DC, se il sum_node ha DC bias)
        # Se i segnali di ingresso sono tutti AC e riferiti a gnd, questo potrebbe non servire.
        self.add_component(Capacitor(name="C_out", node1='sum_node', node2='out', capacitance=1.0e-6))
        self.add_component(Resistor(name="R_out_final_load", node1='out', node2='gnd', resistance=10e3))


# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    # Esempio di utilizzo: due segnali sommati
    passive_mixer = PassiveMixer(num_inputs=2)
    
    print(f"Circuito: {passive_mixer.name}")
    print("Nodi:", passive_mixer.get_node_names())
    print("Componenti:")
    for comp in passive_mixer.components:
        print(f"  - {comp.name}: {type(comp).__name__} ({comp.name})")

    # Aggiungi sorgenti di tensione per testare la somma
    # Nota: Visto che il solver DC non gestisce segnali AC, metteremo tensioni DC per il test DC.
    # In una simulazione transitoria, useresti sorgenti AC.
    passive_mixer.add_component(VoltageSource(name="V_in1", nodes={'pos': 'in1', 'neg': 'gnd'}, voltage=1.0))
    passive_mixer.add_component(VoltageSource(name="V_in2", nodes={'pos': 'in2', 'neg': 'gnd'}, voltage=2.0))
    
    solver = MnaSolver(passive_mixer)

    print("\nSimulazione DC del Passive Mixer...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC ---")
        V_in1 = node_voltages[solver._get_node_index('in1')]
        V_in2 = node_voltages[solver._get_node_index('in2')]
        V_sum_node = node_voltages[solver._get_node_index('sum_node')]
        V_out = node_voltages[solver._get_node_index('out')]

        print(f"V(in1) = {V_in1:.3f} V")
        print(f"V(in2) = {V_in2:.3f} V")
        print(f"V(sum_node) = {V_sum_node:.3f} V")
        print(f"V(out) = {V_out:.3f} V (dovrebbe essere circa V(sum_node) se C_out è ideale in DC)")
        
        # Calcolo atteso per un mixer resistivo
        # V_sum_node = (V_in1 * (1/R_mix1) + V_in2 * (1/R_mix2)) / (1/R_mix1 + 1/R_mix2 + 1/R_sum_load)
        # Se tutte le resistenze sono uguali (R_mix = R_sum_load):
        # V_sum_node = (V_in1 + V_in2) / 3 (per 2 ingressi)
        # V_sum_node = (1.0 + 2.0) / 3 = 1.0 V
        expected_sum = (V_in1 + V_in2) / (passive_mixer.num_inputs + 1) # +1 per R_sum_load
        print(f"Valore atteso di V(sum_node) = {expected_sum:.3f} V")
        if abs(V_sum_node - expected_sum) < 0.01:
            print("Il sommatore passivo funziona come previsto (con attenuazione).")
        else:
            print("ATTENZIONE: Il sommatore passivo non funziona come previsto.")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
