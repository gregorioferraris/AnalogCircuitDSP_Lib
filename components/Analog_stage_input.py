import numpy as np
from circuit_solver.subcircuit import Subcircuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.diode import Diode # Per la distorsione armonica
from components.op_amp import OpAmp # Se vuoi un OpAmp, ma può essere più complesso per MNA non lineare


class AnalogInputStage(Subcircuit):
    """
    Simula uno stadio di ingresso analogico di un convertitore vintage.
    Include:
    - Resistenza di ingresso (Input Impedance)
    - Preamplificazione (con guadagno variabile, implementato con resistori)
    - Una leggera distorsione armonica (usando un diodo in configurazione parallela/serie con resistori)
    - Filtro anti-aliasing passa-basso RC semplice
    """
    def __init__(self, name: str, input_node: str, output_node: str, ground_node: str = '0',
                 input_resistance: float = 10e3,  # Impedenza d'ingresso (10k Ohm)
                 gain_factor: float = 1.0,        # Fattore di guadagno (regolabile da 0 a 1 per atten, >1 per gain)
                 filter_R: float = 1e3,          # Resistenza filtro anti-aliasing
                 filter_C: float = 100e-9,       # Capacità filtro anti-aliasing (100nF)
                 diode_saturation_current: float = 1e-9, # Is per diodo (per distorsione)
                 diode_emission_coefficient: float = 1.0): # N per diodo
        
        super().__init__(name, input_nodes=[input_node], output_nodes=[output_node])
        
        # Nodi interni
        self.gain_node = f"{name}_gain_out"
        self.dist_node = f"{name}_dist_out" # Dopo la distorsione
        self.filter_node = f"{name}_filter_out" # Dopo il filtro

        # Aggiungi i nodi al sottocircuito
        self.add_node(input_node)
        self.add_node(output_node)
        self.add_node(ground_node)
        self.add_node(self.gain_node)
        self.add_node(self.dist_node)
        self.add_node(self.filter_node)

        # 1. Resistenza di Ingresso
        self.add_component(Resistor(f"{name}_Rin", input_node, self.gain_node, input_resistance))

        # 2. Stadio di Guadagno (semplice partitore resistivo per attenuazione/guadagno)
        # Nota: Un guadagno >1.0 richiederà un componente attivo (OpAmp, Transistor)
        # Per semplicità, useremo un partitore per *attenuare* o un semplice resistore.
        # Se gain_factor > 1, questo resistore è solo un passante per l'input_resistance.
        # Per un vero guadagno, useresti un OpAmp.
        if gain_factor < 1.0:
            r_gain_series = input_resistance / gain_factor - input_resistance
            if r_gain_series < 1e-6: r_gain_series = 1e-6 # Evita zero
            self.add_component(Resistor(f"{name}_RGain", self.gain_node, self.dist_node, r_gain_series))
        else: # Se guadagno >= 1, consideralo come passante diretto, gestito da R_in
            self.add_component(Resistor(f"{name}_RPass", self.gain_node, self.dist_node, 1e-6)) # Quasi un corto

        # 3. Distorsione (diodo soft-clipping)
        # Questo è un circuito molto semplice per introdurre non linearità.
        # Il diodo clippa a un certo punto, introducendo armoniche.
        # R_dist_series limita la corrente attraverso il diodo.
        # R_dist_parallel agisce come pull-down per il segnale distorto.
        self.add_component(Resistor(f"{name}_R_dist_series", self.dist_node, f"{name}_d1_in", 1e3))
        self.add_component(Diode(f"{name}_D1", f"{name}_d1_in", ground_node, Is=diode_saturation_current, N=diode_emission_coefficient))
        self.add_component(Resistor(f"{name}_R_dist_parallel", f"{name}_d1_in", self.filter_node, 10e3)) # Connessione al filtro

        # 4. Filtro Anti-Aliasing Passa-Basso (RC)
        self.add_component(Resistor(f"{name}_FilterR", self.filter_node, f"{name}_filter_mid", filter_R))
        self.add_component(Capacitor(f"{name}_FilterC", f"{name}_filter_mid", ground_node, filter_C))
        
        # Output del sottocircuito è il nodo di uscita del filtro
        self.add_component(Resistor(f"{name}_ROut", f"{name}_filter_mid", output_node, 1e-6)) # R quasi 0 per collegare a output_node

        print(f"Stadio di Ingresso Analogico '{name}' creato. Impedenze, guadagno e distorsione modellati.")
