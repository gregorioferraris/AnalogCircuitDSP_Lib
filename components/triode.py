# components/triode.py
import numpy as np
from components.component import Component

class Triode(Component):
    def __init__(self, name: str, plate_node: str, grid_node: str, cathode_node: str,
                 mu: float = 100.0, Kp: float = 1.0, Ex: float = 1.0, Kg1: float = 1.0, Vg_offset: float = 0.0):
        """
        Inizializza un Triodo usando un modello tipo Koren semplificato.
        Args:
            name (str): Nome univoco dell'istanza (es. "V1_preamp").
            plate_node (str): Nome del nodo della placca (anodo).
            grid_node (str): Nome del nodo della griglia.
            cathode_node (str): Nome del nodo del catodo.
            mu (float): Fattore di amplificazione.
            Kp (float): Parametro di permeabilità.
            Ex (float): Esponente.
            Kg1 (float): Parametro di linearità della griglia.
            Vg_offset (float): Offset di tensione per la griglia.
        """
        super().__init__(name, plate_node, grid_node, cathode_node)
        self.pin_names = ('plate', 'grid', 'cathode') # Per mappatura nodi con nome
        self.mu = mu
        self.Kp = Kp
        self.Ex = Ex
        self.Kg1 = Kg1
        self.Vg_offset = Vg_offset # Aggiunto per flessibilità

    def calculate_plate_current(self, Vgk: float, Vpk: float) -> float:
        """
        Calcola la corrente di placca (Ip) per un triodo usando un modello tipo Koren.
        Args:
            Vgk (float): Tensione Griglia-Catodo.
            Vpk (float): Tensione Placca-Catodo.
        Returns:
            float: Corrente di placca in Ampere.
        """
        # Modello Koren semplificato (o una sua variante comune)
        # Assicurati che Vpk sia positivo per evitare logaritmi di numeri negativi
        if Vpk <= 0:
            return 0.0 # Nessuna corrente se la placca non è polarizzata positivamente rispetto al catodo

        # Termine di controllo effettivo della griglia
        # Questo è il cuore del modello Koren per la non linearità
        effective_grid_voltage = Vgk + self.Vg_offset + Vpk / self.mu

        if effective_grid_voltage <= 0:
            return 0.0 # Cutoff

        # Calcolo della corrente di placca
        Ip = self.Kp * (effective_grid_voltage)**self.Ex
        
        # Aggiungi un limite per evitare valori irrealistici di corrente (es. se Kp è molto grande)
        # Questo è un guardrail, non parte del modello fisico puro
        if Ip > 1.0: # Limite arbitrario di 1A per prevenire overflow o instabilità
            Ip = 1.0

        return Ip

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

