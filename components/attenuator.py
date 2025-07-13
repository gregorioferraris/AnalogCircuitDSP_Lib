# components/attenuator.py

from components.resistor import Resistor
from circuit_solver.circuit import Circuit # Importiamo Circuit per coerenza

class Attenuator(Circuit):
    """
    Un attenuatore resistivo variabile simulato come un divisore di tensione fisso.
    Per simularne la variabilità, dovrai modificare r1_val e r2_val esternamente.
    """
    def __init__(self, name, node1, node2, node_tap, attenuation_factor=0.5):
        super().__init__(name)
        
        # attenuation_factor è un valore tra 0 e 1, dove 0 è nessuna uscita e 1 è piena uscita.
        # R1 è la parte superiore del potenziometro, R2 la parte inferiore.
        # Sum_R = R1 + R2
        # V_out = V_in * (R2 / (R1 + R2))
        # attenuation_factor = R2 / (R1 + R2)
        # Per semplicità, useremo una resistenza totale fissa e calcoleremo R1 e R2.
        total_resistance = 100e3 # Ohm (es. un potenziometro da 100k)
        
        # Calcola R2 in base al fattore di attenuazione
        r2_val = total_resistance * attenuation_factor
        r1_val = total_resistance - r2_val

        # Assicurati che le resistenze non siano zero, per evitare problemi nel solver
        if r1_val < 1e-6: r1_val = 1e-6
        if r2_val < 1e-6: r2_val = 1e-6

        self.add_node(node1) # Ingresso dell'attenuatore
        self.add_node(node2) # Massa/riferimento
        self.add_node(node_tap) # Uscita (tap centrale)

        # R1: tra l'ingresso e il tap
        self.add_component(Resistor(name=f"{name}_R1", node1=node1, node2=node_tap, resistance=r1_val))
        # R2: tra il tap e la massa
        self.add_component(Resistor(name=f"{name}_R2", node1=node_tap, node2=node2, resistance=r2_val))
        
        self.input_node = node1
        self.output_node = node_tap
        self.ground_node = node2
        self.r1_val = r1_val
        self.r2_val = r2_val

    def set_attenuation_factor(self, factor):
        """
        Permette di modificare il fattore di attenuazione dinamicamente.
        Questo metodo dovrebbe essere chiamato dal solver durante l'analisi transitoria
        se si vuole simulare un potenziometro che cambia.
        """
        if not (0 <= factor <= 1):
            raise ValueError("Attenuation factor must be between 0 and 1.")
        
        total_resistance = self.components[0].resistance + self.components[1].resistance # Sum of R1 + R2 from init
        
        new_r2_val = total_resistance * factor
        new_r1_val = total_resistance - new_r2_val

        # Evita resistenze troppo piccole
        if new_r1_val < 1e-6: new_r1_val = 1e-6
        if new_r2_val < 1e-6: new_r2_val = 1e-6

        self.components[0].resistance = new_r1_val # R1
        self.components[1].resistance = new_r2_val # R2
        self.r1_val = new_r1_val
        self.r2_val = new_r2_val
        print(f"Attenuator '{self.name}' set to factor {factor:.2f}. R1={new_r1_val:.0f} Ohm, R2={new_r2_val:.0f} Ohm")


# --- Blocco di Test ---
if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    # Esempio di utilizzo: Attenuatore come sottocircuito in un circuito più grande
    class TestCircuit(Circuit):
        def __init__(self):
            super().__init__("Test_Attenuator")
            self.add_node('sig_in')
            self.add_node('attenuated_out')
            
            self.add_component(VoltageSource(name="V_in_test", nodes={'pos': 'sig_in', 'neg': 'gnd'}, voltage=5.0))
            self.attenuator = Attenuator(name="Att_1", node1='sig_in', node2='gnd', node_tap='attenuated_out', attenuation_factor=0.5)
            
            # Aggiungiamo i componenti dell'attenuatore al circuito principale
            for comp in self.attenuator.components:
                self.add_component(comp)

    test_circuit = TestCircuit()
    solver = MnaSolver(test_circuit)
    
    print("\nSimulazione DC con attenuazione 0.5:")
    node_voltages, _ = solver.solve_dc()
    if node_voltages is not None:
        V_in = node_voltages[solver._get_node_index('sig_in')]
        V_out = node_voltages[solver._get_node_index('attenuated_out')]
        print(f"V(sig_in) = {V_in:.3f} V")
        print(f"V(attenuated_out) = {V_out:.3f} V (expected ~ {V_in * 0.5:.3f} V)")

    # Simula il cambio del fattore di attenuazione
    test_circuit.attenuator.set_attenuation_factor(0.25)
    print("\nSimulazione DC con attenuazione 0.25:")
    node_voltages, _ = solver.solve_dc()
    if node_voltages is not None:
        V_in = node_voltages[solver._get_node_index('sig_in')]
        V_out = node_voltages[solver._get_node_index('attenuated_out')]
        print(f"V(sig_in) = {V_in:.3f} V")
        print(f"V(attenuated_out) = {V_out:.3f} V (expected ~ {V_in * 0.25:.3f} V)")
