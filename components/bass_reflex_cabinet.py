# components/bass_reflex_cabinet.py
import numpy as np
from components.component import Component

class BassReflexCabinet(Component):
    def __init__(self, name: str, mech_velocity_input_node: str,
                 box_volume: float, port_area: float, port_length: float,
                 air_density: float = 1.225, speed_of_sound: float = 343.0, port_loss_factor: float = 0.5):
        """
        Inizializza un modello di cabinet bass reflex.
        Include la cedevolezza della cassa e la risonanza del condotto (port).
        Args:
            name (str): Nome univoco dell'istanza (es. "BR1").
            mech_velocity_input_node (str): Nodo che si connette alla velocità meccanica del cono (tensione analogica).
            box_volume (float): Volume interno della cassa in metri cubi (m^3).
            port_area (float): Area della sezione trasversale del condotto (m^2).
            port_length (float): Lunghezza del condotto (m).
            air_density (float): Densità dell'aria (kg/m^3).
            speed_of_sound (float): Velocità del suono nell'aria (m/s).
            port_loss_factor (float): Fattore per modellare le perdite nel condotto (adimensionale).
        """
        super().__init__(name, mech_velocity_input_node)
        self.pin_names = ('mech_velocity_in',) # Nomi dei pin per mappatura

        self.box_volume = box_volume
        self.port_area = port_area
        self.port_length = port_length
        self.air_density = air_density
        self.speed_of_sound = speed_of_sound
        self.port_loss_factor = port_loss_factor

        # Calcola la cedevolezza acustica della cassa (analogia elettrica: Capacità)
        self.C_box = self.box_volume / (self.air_density * self.speed_of_sound**2)

        # Calcola la massa acustica del condotto (analogia elettrica: Induttanza)
        # Aggiustamento per la lunghezza effettiva del port (end correction)
        effective_port_length = self.port_length + 0.85 * np.sqrt(self.port_area / np.pi) # Approx. per un'estremità libera
        self.L_port = self.air_density * effective_port_length / self.port_area

        # Calcola la resistenza acustica del condotto (analogia elettrica: Resistenza)
        # Questo è un modello semplificato per le perdite del port.
        self.R_port = self.port_loss_factor * (self.air_density * self.speed_of_sound / self.port_area)

        # Variabili di stato per C_box, L_port
        self.v_Cbox_prev = 0.0
        self.i_Cbox_prev = 0.0
        self.v_Lport_prev = 0.0
        self.i_Lport_prev = 0.0

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi del cabinet bass reflex alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        Il modello prevede:
        - C_box (Capacitor) tra mech_velocity_input_node e ground.
        - L_port (Inductor) in serie con R_port (Resistor) tra mech_velocity_input_node e ground.
          Questo crea un ramo risonante.

        Questo è un placeholder concettuale. L'implementazione completa richiederebbe:
        - Creazione di nodi MNA interni per il ramo del port (es. nodo tra L_port e R_port).
        - Utilizzo dei modelli di Resistor, Capacitor, Inductor per i rispettivi elementi.
        - Connessione di questi elementi ai nodi appropriati (mech_velocity_input_node, nodi interni, ground).
        - Gestione dello stato interno (v_prev, i_prev) per i componenti dinamici equivalenti.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        mech_velocity_in_id = self.node_ids['mech_velocity_in']
        ground_id = 0 # Il ground acustico è il ground elettrico

        # --- Contributi di C_box (Capacità della cassa) ---
        G_eq_Cbox = 2.0 * self.C_box / dt
        if mech_velocity_in_id != 0: stamp_A[mech_velocity_in_id, mech_velocity_in_id] += G_eq_Cbox
        V_Cbox_prev = prev_solution[mech_velocity_in_id] - prev_solution[ground_id]
        i_eq_Cbox = G_eq_Cbox * V_Cbox_prev + self.i_Cbox_prev
        if mech_velocity_in_id != 0: stamp_B[mech_velocity_in_id] -= i_eq_Cbox

        # --- Contributi del ramo del port (L_port in serie con R_port) ---
        # Questo richiede un nodo interno per modellare correttamente la serie L-R in parallelo con C_box.
        # Per semplicità in get_stamps, lo tratteremo come un'impedenza equivalente.
        # In una simulazione MNA completa, dovresti aggiungere un nodo interno e gli stamps per L e R separatamente.
        
        # Per ora, come placeholder, non aggiungiamo stamps diretti per il ramo L_port/R_port qui,
        # poiché la loro interazione completa richiede nodi interni che non sono gestiti da questa get_stamps.
        # La loro gestione nello _system_equations del MnaSolver o tramite un sottocircuito sarebbe più adatta.

        return stamp_A, stamp_B

    def update_state(self, v_Cbox_curr: float, i_Cbox_curr: float, v_Lport_curr: float, i_Lport_curr: float):
        """
        Aggiorna lo stato interno dei componenti dinamici equivalenti.
        Questo metodo dovrebbe essere chiamato dal MnaSolver per i componenti dinamici.
        """
        self.v_Cbox_prev = v_Cbox_curr
        self.i_Cbox_prev = i_Cbox_curr
        self.v_Lport_prev = v_Lport_curr
        self.i_Lport_prev = i_Lport_curr

    def get_acoustic_impedance_analogy(self, frequency: float) -> complex:
        """
        Calcola l'impedenza acustica equivalente (in analogia elettrica) del cabinet bass reflex.
        Args:
            frequency (float): Frequenza in Hz.
        Returns:
            complex: Impedenza acustica in analogia elettrica.
        """
        omega = 2 * np.pi * frequency
        
        # Impedenza della capacità della cassa (C_box)
        Z_Cbox = 1.0 / (1j * omega * self.C_box) if omega != 0 else np.inf
        
        # Impedenza dell'induttanza del port (L_port)
        Z_Lport = 1j * omega * self.L_port
        
        # Impedenza della resistenza del port (R_port)
        Z_Rport = self.R_port
        
        # Il ramo del port è L_port in serie con R_port
        Z_port_branch = Z_Lport + Z_Rport
        
        # Il cabinet bass reflex è C_box in parallelo con il ramo del port
        # 1/Z_total = 1/Z_Cbox + 1/Z_port_branch
        Y_total = 0.0
        if Z_Cbox != 0:
            Y_total += 1.0 / Z_Cbox
        if Z_port_branch != 0:
            Y_total += 1.0 / Z_port_branch
        
        return 1.0 / Y_total if Y_total != 0 else np.inf

