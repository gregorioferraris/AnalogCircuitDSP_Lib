# components/closed_box_cabinet.py
import numpy as np
from components.component import Component

class ClosedBoxCabinet(Component):
    def __init__(self, name: str, mech_velocity_input_node: str,
                 box_volume: float, air_density: float = 1.225, speed_of_sound: float = 343.0):
        """
        Inizializza un modello di cabinet a cassa chiusa.
        Modella la cedevolezza acustica dell'aria all'interno della cassa.
        Args:
            name (str): Nome univoco dell'istanza (es. "CB1").
            mech_velocity_input_node (str): Nodo che si connette alla velocità meccanica del cono (tensione analogica).
            box_volume (float): Volume interno della cassa in metri cubi (m^3).
            air_density (float): Densità dell'aria (kg/m^3).
            speed_of_sound (float): Velocità del suono nell'aria (m/s).
        """
        super().__init__(name, mech_velocity_input_node)
        self.pin_names = ('mech_velocity_in',) # Nomi dei pin per mappatura

        self.box_volume = box_volume
        self.air_density = air_density
        self.speed_of_sound = speed_of_sound

        # Calcola la cedevolezza acustica della cassa (analogia elettrica: Capacità)
        # C_box = Vb / (rho * c^2)
        self.C_box = self.box_volume / (self.air_density * self.speed_of_sound**2)

        # Variabili di stato per C_box
        self.v_Cbox_prev = 0.0
        self.i_Cbox_prev = 0.0

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi del cabinet a cassa chiusa alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        Il cabinet agisce come una capacità acustica (C_box) connessa tra il nodo di velocità meccanica
        e il riferimento (ground acustico, che è lo stesso ground elettrico in questa analogia).
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        mech_velocity_in_id = self.node_ids['mech_velocity_in']
        ground_id = 0 # Il ground acustico è il ground elettrico

        # Contributi di C_box
        G_eq_Cbox = 2.0 * self.C_box / dt
        
        if mech_velocity_in_id != 0: stamp_A[mech_velocity_in_id, mech_velocity_in_id] += G_eq_Cbox
        # Nota: il ground (ID 0) non ha una riga/colonna esplicita nella matrice MNA risolta,
        # ma è il riferimento per le tensioni.

        V_Cbox_prev = prev_solution[mech_velocity_in_id] - prev_solution[ground_id]
        i_eq_Cbox = G_eq_Cbox * V_Cbox_prev + self.i_Cbox_prev
        
        if mech_velocity_in_id != 0: stamp_B[mech_velocity_in_id] -= i_eq_Cbox

        return stamp_A, stamp_B

    def update_state(self, v_curr: float, i_curr: float):
        """
        Aggiorna lo stato interno della capacità acustica per il prossimo passo temporale.
        Questo metodo dovrebbe essere chiamato dal MnaSolver per i componenti dinamici.
        """
        self.v_Cbox_prev = v_curr
        self.i_Cbox_prev = i_curr

    def get_acoustic_impedance_analogy(self, frequency: float) -> complex:
        """
        Calcola l'impedenza acustica equivalente (in analogia elettrica) del cabinet a cassa chiusa.
        Args:
            frequency (float): Frequenza in Hz.
        Returns:
            complex: Impedenza acustica in analogia elettrica.
        """
        omega = 2 * np.pi * frequency
        # Un condensatore in analogia elettrica ha impedenza 1/(jwC)
        if omega == 0:
            return np.inf # Circuito aperto a DC
        return 1.0 / (1j * omega * self.C_box)

