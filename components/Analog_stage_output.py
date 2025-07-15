import numpy as np
from circuit_solver.subcircuit import Subcircuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.diode import Diode # Per la distorsione/limitazione


class AnalogOutputStage(Subcircuit):
    """
    Simula uno stadio di uscita analogico di un convertitore vintage.
    Include:
    - Filtro di ricostruzione passa-basso RC
    - Resistenza di uscita (Output Impedance)
    - Possibile limitazione/saturazione (via diodi)
    """
    def __init__(self, name: str, input_node: str, output_node: str, ground_node: str = '0',
                 filter_R: float = 1e3,          # Resistenza filtro di ricostruzione
                 filter_C: float = 100e-9,       # Capacità filtro di ricostruzione
                 output_resistance: float = 50.0, # Impedenza d'uscita (50 Ohm)
                 clipper_threshold: float = 0.7, # Soglia approssimativa di clipping per diodi
                 diode_saturation_current: float = 1e-9, # Is per diodo clipper
                 diode_emission_coefficient: float = 1.0): # N per diodo clipper
        
        super().__init__(name, input_nodes=[input_node], output_nodes=[output_node])

        # Nodi interni
        self.filter_in_node = f"{name}_filter_in"
        self.clipper_node = f"{name}_clipper_out"
        
        # Aggiungi i nodi al sottocircuito
        self.add_node(input_node)
        self.add_node(output_node)
        self.add_node(ground_node)
        self.add_node(self.filter_in_node)
        self.add_node(self.clipper_node)

        # 1. Filtro di Ricostruzione Passa-Basso (RC)
        # Questo serve a "lisciare" il segnale dopo il DAC digitale (non simulato qui)
        self.add_component(Resistor(f"{name}_FilterR_in", input_node, self.filter_in_node, filter_R))
        self.add_component(Capacitor(f"{name}_FilterC_out", self.filter_in_node, ground_node, filter_C))

        # 2. Stadio di Clipper/Saturazione (diode limiter)
        # Per simulare una leggera compressione o saturazione.
        # Nota: Questo è un clipper duro. Potresti voler usare qualcosa di più soft
        # come un OpAmp con feedback non lineare per un effetto più graduale.
        # Qui usiamo un diodo in parallelo per tagliare il segnale a un certo livello.
        # D2 è un diodo che clippa il segnale positivo.
        # D3 è un diodo che clippa il segnale negativo (collegato al contrario).
        # R_clipper_series limita la corrente attraverso i diodi.
        self.add_component(Resistor(f"{name}_R_clipper_series", self.filter_in_node, self.clipper_node, 100)) # Limita corrente
        
        # Diodi per clipping bidirezionale
        # Per controllare la soglia di clip:
        # Usa un valore alto di Is per simulare una tensione di soglia più bassa (vicino a 0.2-0.3V)
        # Oppure metti più diodi in serie per aumentare la soglia di ~0.7V per diodo.
        # Per clipper_threshold: useremo Rs per approssimare il livello di clipping desiderato.
        # Un diodo singolo si attiva intorno a 0.7V. Per clip a 0.3V, Is dev'essere più alta (es. 1e-6).
        # Se vogliamo una soglia specifica, dovremmo modellare con una sorgente di tensione in serie ai diodi.
        
        # Per simulare una soglia di clipping approssimativa con un diodo singolo (o parallelo)
        # usiamo Is per "spostare" la curva I-V. Valori più alti di Is significano che il diodo
        # conduce di più a tensioni più basse.
        
        # Diodo per clipping positivo (taglia i picchi sopra la soglia del diodo)
        self.add_component(Diode(f"{name}_D_pos", self.clipper_node, ground_node,
                                  Is=diode_saturation_current * (10**(clipper_threshold/0.026)) , # Ajusta Is per soglia desiderata (V_T=26mV)
                                  N=diode_emission_coefficient))
        # Diodo per clipping negativo (taglia i picchi sotto la soglia negativa del diodo)
        self.add_component(Diode(f"{name}_D_neg", ground_node, self.clipper_node,
                                  Is=diode_saturation_current * (10**(clipper_threshold/0.026)),
                                  N=diode_emission_coefficient))

        # 3. Resistenza di Uscita
        self.add_component(Resistor(f"{name}_Rout", self.clipper_node, output_node, output_resistance))

        print(f"Stadio di Uscita Analogico '{name}' creato. Filtro, impedenza d'uscita e clipping modellati.")
