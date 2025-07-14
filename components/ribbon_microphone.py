# components/ribbon_microphone.py
import numpy as np
from components.component import Component

class RibbonMicrophone(Component):
    def __init__(self, name: str, output_node: str, ground_node: str,
                 ribbon_area: float = 5e-5, ribbon_mass: float = 1e-6, ribbon_compliance: float = 1e-5,
                 mechanical_resistance: float = 0.01, Bl: float = 0.5, transformer_ratio: float = 30.0, output_impedance: float = 300.0):
        """
        Inizializza un modello semplificato di microfono a nastro.
        Args:
            name (str): Nome univoco dell'istanza (es. "RMIC1").
            output_node (str): Nodo di uscita del segnale elettrico.
            ground_node (str): Nodo di riferimento (ground).
            ribbon_area (float): Area effettiva del nastro (m^2).
            ribbon_mass (float): Massa del nastro (kg).
            ribbon_compliance (float): Cedevolezza meccanica del nastro (m/N).
            mechanical_resistance (float): Resistenza meccanica delle perdite (Ns/m).
            Bl (float): Fattore di forza (Tesla-metri o N/A).
            transformer_ratio (float): Rapporto di spire del trasformatore di uscita (primario:secondario).
            output_impedance (float): Impedenza di uscita nominale del microfono (Ohm).
        """
        super().__init__(name, output_node, ground_node)
        self.pin_names = ('output', 'ground')

        self.ribbon_area = ribbon_area
        self.ribbon_mass = ribbon_mass
        self.ribbon_compliance = ribbon_compliance
        self.mechanical_resistance = mechanical_resistance
        self.Bl = Bl
        self.transformer_ratio = transformer_ratio # N_sec / N_prim
        self.output_impedance = output_impedance # Impedenza vista dall'esterno

        # Parametri meccanici equivalenti (analogia elettrica)
        self.L_mech_eq = self.ribbon_mass
        self.C_mech_eq = self.ribbon_compliance
        self.R_mech_eq = self.mechanical_resistance

    def calculate_output_voltage(self, acoustic_pressure: float) -> float:
        """
        Calcola la tensione di uscita del microfono a nastro data la pressione acustica.
        Modello molto semplificato.
        Args:
            acoustic_pressure (float): Pressione sonora in Pascal (Pa).
        Returns:
            float: Tensione di uscita in Volt.
        """
        # Forza sul nastro: F = P * A
        force = acoustic_pressure * self.ribbon_area

        # La tensione generata dal nastro è proporzionale alla sua velocità (V_ribbon = Bl * v_ribbon)
        # Questa tensione è poi trasformata dal trasformatore di uscita.
        # Sensibilità tipica: 1-2 mV/Pa (0.001 - 0.002 V/Pa)
        sensitivity_factor = (self.Bl * self.ribbon_area / (self.mechanical_resistance + 1e-6)) * self.transformer_ratio # Semplificazione

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
        Calcola la risposta in frequenza del microfono a nastro.
        Include la risonanza meccanica del nastro e l'effetto del trasformatore.
        Args:
            frequencies (np.ndarray): Frequenze in Hz.
        Returns:
            np.ndarray: Risposta in frequenza complessa (V/Pa).
        """
        omega = 2 * np.pi * frequencies
        
        # Impedenza meccanica del nastro
        Z_mech = self.R_mech_eq + 1j * omega * self.L_mech_eq + 1.0 / (1j * omega * self.C_mech_eq)
        
        # La tensione generata dal nastro è V_ribbon = Bl * v_ribbon
        # dove v_ribbon = Forza / Z_mech = (P_ac * Area_ribbon) / Z_mech
        # Quindi V_ribbon = (Bl * P_ac * Area_ribbon) / Z_mech
        
        # Il trasformatore di uscita ha un rapporto N_sec/N_prim = transformer_ratio
        # La tensione di uscita è V_out = V_ribbon * transformer_ratio
        
        # Risposta in frequenza: V_out / P_ac
        response = (self.Bl * self.ribbon_area * self.transformer_ratio) / Z_mech if Z_mech != 0 else np.inf
        
        # I microfoni a nastro hanno tipicamente un roll-off in alta frequenza
        # dovuto all'induttanza del trasformatore o alla massa del nastro.
        # Questo modello semplificato non lo cattura completamente.
        
        return response

