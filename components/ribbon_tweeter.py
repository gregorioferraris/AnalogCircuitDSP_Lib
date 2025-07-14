# components/ribbon_tweeter.py
import numpy as np
from components.component import Component

class RibbonTweeter(Component):
    def __init__(self, name: str, elec_plus_node: str, elec_minus_node: str, mech_velocity_node: str,
                 Re: float, Le: float, Bl: float, Mms: float, Cms: float, Rms: float, Sd: float, transformer_ratio: float):
        """
        Inizializza un modello di tweeter a nastro basato su parametri Thiele-Small.
        Include un trasformatore di accoppiamento per l'impedenza.
        Args:
            name (str): Nome univoco dell'istanza (es. "RTW1").
            elec_plus_node (str): Nodo elettrico positivo del trasformatore.
            elec_minus_node (str): Nodo elettrico negativo del trasformatore.
            mech_velocity_node (str): Nodo che rappresenta la velocità meccanica del nastro.
            Re (float): Resistenza DC del nastro (molto bassa, Ohm).
            Le (float): Induttanza del nastro (molto bassa, Henry).
            Bl (float): Fattore di forza (molto basso, Tesla-metri o N/A).
            Mms (float): Massa meccanica in movimento del nastro (molto bassa, kg).
            Cms (float): Cedevolezza meccanica della sospensione del nastro (molto alta, m/N).
            Rms (float): Resistenza meccanica delle perdite (Ns/m).
            Sd (float): Area effettiva del nastro (m^2).
            transformer_ratio (float): Rapporto di spire del trasformatore di accoppiamento (primario:secondario).
        """
        super().__init__(name, elec_plus_node, elec_minus_node, mech_velocity_node)
        self.pin_names = ('elec_plus', 'elec_minus', 'mech_velocity_out')

        # Parametri elettrici del nastro (lato primario del trasformatore)
        self.Re_ribbon = Re
        self.Le_ribbon = Le

        # Parametri meccanici (analogia elettrica)
        self.L_mech = Mms
        self.C_mech = Cms
        self.R_mech = Rms
        self.Sd = Sd

        # Fattore di accoppiamento elettro-meccanico del nastro
        self.Bl_ribbon = Bl

        # Trasformatore di accoppiamento
        self.transformer_ratio = transformer_ratio # N_prim / N_sec (per impedenza)

        # Variabili di stato per i componenti dinamici (se gestite internamente)
        self.v_Le_prev = 0.0
        self.i_Le_prev = 0.0
        self.v_Lmech_prev = 0.0
        self.i_Lmech_prev = 0.0
        self.v_Cmech_prev = 0.0
        self.i_Cmech_prev = 0.0

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi del tweeter a nastro alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        Questo metodo implementa un modello a elementi concentrati.
        Simile a SpeakerDriver, ma con il trasformatore.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        elec_plus_id = self.node_ids['elec_plus']
        elec_minus_id = self.node_ids['elec_minus']
        mech_velocity_out_id = self.node_ids['mech_velocity_out']
        ground_id = 0

        # --- Parte Elettrica (vista dal lato primario del trasformatore) ---
        # L'impedenza del nastro è riflessa sul primario del trasformatore.
        # Z_ribbon_elec = Re_ribbon + jwLe_ribbon
        # Z_ribbon_mech_reflected = (Bl_ribbon^2) / Z_mech_ribbon
        # Z_total_ribbon_side = Z_ribbon_elec + Z_ribbon_mech_reflected
        
        # Questa impedenza è poi trasformata dal trasformatore al lato di ingresso.
        # Z_input = Z_total_ribbon_side * (transformer_ratio)^2
        
        # Per MNA, questo è complesso da implementare direttamente in get_stamps.
        # Richiederebbe nodi interni per il trasformatore e per il lato nastro.
        # Per ora, è un placeholder.

        return stamp_A, stamp_B

    def update_state(self, v_Le_curr: float, i_Le_curr: float, v_Lmech_curr: float, i_Lmech_curr: float, v_Cmech_curr: float, i_Cmech_curr: float):
        """
        Aggiorna lo stato interno dei componenti dinamici equivalenti.
        """
        self.v_Le_prev = v_Le_curr
        self.i_Le_prev = i_Le_curr
        self.v_Lmech_prev = v_Lmech_curr
        self.i_Lmech_prev = i_Lmech_curr
        self.v_Cmech_prev = v_Cmech_curr
        self.i_Cmech_prev = i_Cmech_curr

    def get_electrical_impedance(self, frequency: float) -> complex:
        """
        Calcola l'impedenza elettrica vista dai terminali di ingresso del tweeter (lato primario del trasformatore).
        Args:
            frequency (float): Frequenza in Hz.
        Returns:
            complex: Impedenza elettrica.
        """
        omega = 2 * np.pi * frequency
        
        # Impedenza elettrica del nastro
        Z_ribbon_elec = self.Re_ribbon + 1j * omega * self.Le_ribbon
        
        # Impedenza meccanica del nastro (analogia elettrica)
        Z_mech_ribbon = self.R_mech + 1j * omega * self.L_mech + 1.0 / (1j * omega * self.C_mech)
        
        # Impedenza elettrica riflessa dalla parte meccanica sul lato nastro
        Z_reflected_ribbon_side = (self.Bl_ribbon**2) / Z_mech_ribbon if Z_mech_ribbon != 0 else np.inf
        
        # Impedenza totale sul lato nastro del trasformatore
        Z_total_ribbon_side = Z_ribbon_elec + Z_reflected_ribbon_side
        
        # Trasforma l'impedenza al lato primario del trasformatore
        # Z_primary = Z_secondary * (N_prim / N_sec)^2
        # Qui transformer_ratio = N_prim / N_sec
        Z_input = Z_total_ribbon_side * (self.transformer_ratio**2)
        
        return Z_input

    def get_acoustic_output_gain(self, frequency: float, input_voltage: float = 1.0) -> complex:
        """
        Calcola il guadagno di uscita acustica (pressione sonora per tensione in ingresso) per il tweeter a nastro.
        Args:
            frequency (float): Frequenza in Hz.
            input_voltage (float): Tensione RMS di ingresso al tweeter.
        Returns:
            complex: Guadagno acustico (Pa/V).
        """
        omega = 2 * np.pi * frequency
        
        # Impedenza elettrica totale vista dall'ingresso (già calcolata)
        Z_input = self.get_electrical_impedance(frequency)
        
        # Corrente sul lato primario del trasformatore: I_prim = input_voltage / Z_input
        # Corrente sul lato nastro: I_ribbon = I_prim * transformer_ratio
        I_ribbon = (input_voltage / Z_input) * self.transformer_ratio if Z_input != 0 else 0.0
        
        # Forza sul nastro: F_ribbon = Bl_ribbon * I_ribbon
        F_ribbon = self.Bl_ribbon * I_ribbon
        
        # Velocità del nastro: v_ribbon = F_ribbon / Z_mech_ribbon
        Z_mech_ribbon = self.R_mech + 1j * omega * self.L_mech + 1.0 / (1j * omega * self.C_mech)
        v_ribbon = F_ribbon / Z_mech_ribbon if Z_mech_ribbon != 0 else 0.0
        
        # Pressione sonora (proporzionale a velocità * area)
        # P_ac = rho * c * v_ribbon * Sd (semplificato)
        # Qui, restituiamo la velocità del nastro normalizzata per l'input voltage e l'area del nastro.
        # Questo è il guadagno dalla tensione elettrica alla velocità del nastro.
        acoustic_output_gain = (v_ribbon * self.Sd) / input_voltage if input_voltage != 0 else 0.0
        
        return acoustic_output_gain

