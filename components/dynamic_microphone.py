# components/dynamic_microphone.py
import numpy as np
from components.component import Component

class DynamicMicrophone(Component):
    def __init__(self, name: str, output_node: str, ground_node: str,
                 diaphragm_area: float = 1e-4, coil_mass: float = 1e-3, suspension_compliance: float = 1e-6,
                 mechanical_resistance: float = 0.5, Bl: float = 5.0, coil_resistance: float = 150.0, coil_inductance: float = 0.5e-3):
        """
        Inizializza un modello semplificato di microfono dinamico (a bobina mobile).
        Args:
            name (str): Nome univoco dell'istanza (es. "DMIC1").
            output_node (str): Nodo di uscita del segnale elettrico.
            ground_node (str): Nodo di riferimento (ground).
            diaphragm_area (float): Area effettiva del diaframma (m^2).
            coil_mass (float): Massa della bobina mobile e diaframma (kg).
            suspension_compliance (float): Cedevolezza meccanica della sospensione (m/N).
            mechanical_resistance (float): Resistenza meccanica delle perdite (Ns/m).
            Bl (float): Fattore di forza (Tesla-metri o N/A).
            coil_resistance (float): Resistenza DC della bobina mobile (Ohm).
            coil_inductance (float): Induttanza della bobina mobile (Henry).
        """
        super().__init__(name, output_node, ground_node)
        self.pin_names = ('output', 'ground')

        self.diaphragm_area = diaphragm_area
        self.coil_mass = coil_mass
        self.suspension_compliance = suspension_compliance
        self.mechanical_resistance = mechanical_resistance
        self.Bl = Bl
        self.coil_resistance = coil_resistance
        self.coil_inductance = coil_inductance

        # Parametri meccanici equivalenti (analogia elettrica)
        self.L_mech_eq = self.coil_mass
        self.C_mech_eq = self.suspension_compliance
        self.R_mech_eq = self.mechanical_resistance

    def calculate_output_voltage(self, acoustic_pressure: float) -> float:
        """
        Calcola la tensione di uscita del microfono dinamico data la pressione acustica.
        Questo è un modello molto semplificato.
        Args:
            acoustic_pressure (float): Pressione sonora in Pascal (Pa).
        Returns:
            float: Tensione di uscita in Volt.
        """
        # Forza sul diaframma: F = P * A
        force = acoustic_pressure * self.diaphragm_area

        # La tensione di uscita è proporzionale alla velocità della bobina (V = Bl * v)
        # e la velocità dipende dalla forza e dall'impedenza meccanica.
        # Per ora, un modello lineare basato sulla sensibilità.
        # Sensibilità tipica: 1-3 mV/Pa (0.001 - 0.003 V/Pa)
        sensitivity_factor = self.Bl * self.diaphragm_area / (self.mechanical_resistance + 1e-6) # Semplificazione

        output_voltage = acoustic_pressure * sensitivity_factor
        
        return output_voltage

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i microfoni, get_stamps non contribuisce direttamente alla matrice MNA
        nel modo tradizionale, poiché sono generatori di segnale basati su un input acustico.
        La loro uscita è una sorgente di tensione (o corrente) controllata esternamente.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

    def get_frequency_response(self, frequencies: np.ndarray) -> np.ndarray:
        """
        Calcola la risposta in frequenza del microfono dinamico.
        Include la risonanza meccanica e l'impedenza elettrica della bobina.
        Args:
            frequencies (np.ndarray): Frequenze in Hz.
        Returns:
            np.ndarray: Risposta in frequenza complessa (V/Pa).
        """
        omega = 2 * np.pi * frequencies
        
        # Impedenza meccanica del sistema (bobina+diaframma+sospensione)
        Z_mech = self.R_mech_eq + 1j * omega * self.L_mech_eq + 1.0 / (1j * omega * self.C_mech_eq)
        
        # Impedenza elettrica della bobina
        Z_elec_coil = self.coil_resistance + 1j * omega * self.coil_inductance
        
        # Funzione di trasferimento dalla pressione acustica alla tensione di uscita
        # V_out / P_ac = (Bl * Area_diaframma) / (Z_mech + (Bl^2 / Z_elec_coil))
        
        # Termine di impedenza meccanica riflessa dalla parte elettrica
        Z_reflected_mech = (self.Bl**2) / Z_elec_coil if Z_elec_coil != 0 else np.inf
        
        denominator = Z_mech + Z_reflected_mech
        
        # Sensibilità (V/Pa)
        response = (self.Bl * self.diaphragm_area) / denominator if denominator != 0 else np.inf
        
        return response

