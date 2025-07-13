# circuits/TubePreampCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
#from components.triode import Triode # Assumi che il tuo Triode sia qui e abbia metodi simili
                                      # a calculate_plate_current e calculate_jacobian_elements
# Per ora uso un placeholder, dovrai importare il tuo Triode reale
# Se il tuo Triodo è molto diverso, dovremo adattare.
# Assumo che abbia un __init__(name, mu, K, Vg_cutoff, rp) e set_nodes(anode, grid, cathode)
class TriodePlaceholder:
    def __init__(self, name="Triode", mu=100.0, K=100.0e-6, Vg_cutoff=-1.5, rp=62500.0):
        self.name = name
        self.mu = float(mu)
        self.K = float(K)
        self.Vg_cutoff = float(Vg_cutoff)
        self.rp = float(rp)
        self.nodes = {}

    def set_nodes(self, anode_node_id, grid_node_id, cathode_node_id):
        self.nodes['anode'] = anode_node_id
        self.nodes['grid'] = grid_node_id
        self.nodes['cathode'] = cathode_node_id

    def calculate_plate_current(self, v_gk, v_pk):
        V_eff = v_gk + (v_pk / self.mu)
        if V_eff <= self.Vg_cutoff: return 0.0
        return max(0.0, self.K * (V_eff - self.Vg_cutoff)**1.5)

    def calculate_jacobian_elements(self, v_anode, v_grid, v_cathode):
        v_gk = v_grid - v_cathode
        v_pk = v_anode - v_cathode
        gm = numerical_jacobian(lambda v: self.calculate_plate_current(v, v_pk), v_gk)
        gp = numerical_jacobian(lambda v: self.calculate_plate_current(v_gk, v), v_pk)
        jacobian = np.zeros((3, 3))
        jacobian[0, 0] = gp; jacobian[0, 1] = gm; jacobian[0, 2] = -(gm + gp)
        jacobian[1, :] = 0.0
        jacobian[2, 0] = -gp; jacobian[2, 1] = -gm; jacobian[2, 2] = (gm + gp)
        return jacobian

    def __str__(self): return f"Triode({self.name})"

# Sostituisci "TriodePlaceholder" con "Triode" una volta che avrai il tuo vero modulo Triode.
# from components.triode import Triode # DECOMMENTA QUESTA RIGA E RIMUOVI IL PLACEHOLDER SOPRA!

class TubePreampCircuit(Circuit):
    """
    Circuito di un semplice preamplificatore audio basato su una valvola Triodo.
    Configurazione a catodo comune (Common Cathode) con bias resistivo.
    """
    def __init__(self, name="Tube_Preamp",
                 tube_params=None, # Parametri del Triodo
                 R_grid=1.0e6,      # Resistenza di griglia (Grid Leak Resistor)
                 R_cathode=1500.0,  # Resistenza di catodo (Cathode Resistor for bias)
                 C_cathode=22.0e-6, # Condensatore di bypass del catodo
                 R_anode=100000.0,  # Resistenza di placca (Anode Load Resistor)
                 C_coupling_in=0.1e-6, # Condensatore di accoppiamento in ingresso
                 C_coupling_out=0.1e-6, # Condensatore di accoppiamento in uscita
                 V_power_supply=250.0, # Tensione di alimentazione DC (B+)
                 sample_rate=48000):
        super().__init__(name)
        self.sample_rate = sample_rate

        self.tube_params = tube_params if tube_params is not None else {}
        self.R_grid_val = float(R_grid)
        self.R_cathode_val = float(R_cathode)
        self.C_cathode_val = float(C_cathode)
        self.R_anode_val = float(R_anode)
        self.C_coupling_in_val = float(C_coupling_in)
        self.C_coupling_out_val = float(C_coupling_out)
        self.V_power_supply = float(V_power_supply)

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per il preamplificatore."""
        self.add_node("Input")
        self.add_node("Output")
        self.add_node("V_Power_Supply") # Nodo per l'alimentazione B+

        # Nodi della Valvola
        self.add_node("Tube_Anode")
        self.add_node("Tube_Grid")
        self.add_node("Tube_Cathode")

        # Nodi intermedi
        self.add_node("Input_Coupling_Node") # Tra C_coupling_in e R_grid
        self.add_node("Output_Coupling_Node") # Tra Anode e C_coupling_out

    def _add_components(self):
        """Aggiunge i componenti (Valvola, R, C)."""
        # Sorgente di alimentazione B+
        self.add_voltage_source(self.V_power_supply, "V_Power_Supply", "GND", name="B_Plus_Supply")

        # Valvola Triodo
        self.tube = TriodePlaceholder(name="Preamp_Tube", **self.tube_params) # Usa TriodePlaceholder per ora
        self.tube.set_nodes(anode_node_id=self.get_node_id("Tube_Anode"),
                            grid_node_id=self.get_node_id("Tube_Grid"),
                            cathode_node_id=self.get_node_id("Tube_Cathode"))
        self.nonlinear_components.append(self.tube) # Aggiungi alla lista dei non lineari

        # Componenti di ingresso
        self.add_component(Capacitor(self.C_coupling_in_val, sample_rate=self.sample_rate),
                           "Input", "Input_Coupling_Node", name="C_coupling_in")
        self.add_component(Resistor(self.R_grid_val),
                           "Input_Coupling_Node", "GND", name="R_grid")
        
        # Connessione della griglia della valvola
        self.connect_nodes("Tube_Grid", "Input_Coupling_Node")

        # Componenti del catodo (bias e bypass)
        self.add_component(Resistor(self.R_cathode_val),
                           "Tube_Cathode", "GND", name="R_cathode")
        self.add_component(Capacitor(self.C_cathode_val, sample_rate=self.sample_rate),
                           "Tube_Cathode", "GND", name="C_cathode") # In parallelo a R_cathode

        # Componenti dell'anodo (carico e accoppiamento uscita)
        self.add_component(Resistor(self.R_anode_val),
                           "V_Power_Supply", "Tube_Anode", name="R_anode")
        self.add_component(Capacitor(self.C_coupling_out_val, sample_rate=self.sample_rate),
                           "Tube_Anode", "Output_Coupling_Node", name="C_coupling_out")

    def _connect_nodes(self):
        """Connette i nodi come richiesto."""
        self.connect_nodes("Output", "Output_Coupling_Node")

    def get_input_node(self):
        return self.get_node_id("Input")

    def get_output_node(self):
        return self.get_node_id("Output")
