# components/zener_diode.py
import numpy as np
from components.component import Component

class ZenerDiode(Component):
    def __init__(self, name: str, anode_node: str, cathode_node: str,
                 Is: float = 1e-14, N: float = 1.0, Vt: float = 0.0258, Vz: float = 5.1, Iz: float = 1e-3, Rz: float = 10.0):
        """
        Inizializza un diodo Zener.
        Modella la conduzione in polarizzazione diretta (come un diodo normale)
        e la rottura in polarizzazione inversa alla tensione Zener.
        Args:
            name (str): Nome univoco dell'istanza (es. "ZD1").
            anode_node (str): Nome del nodo dell'anodo.
            cathode_node (str): Nome del nodo del catodo.
            Is (float): Corrente di saturazione inversa (per la parte diodo normale).
            N (float): Fattore di idealità (per la parte diodo normale).
            Vt (float): Tensione termica (kT/q).
            Vz (float): Tensione Zener (tensione di rottura inversa).
            Iz (float): Corrente di test Zener (corrente alla quale Vz è specificata).
            Rz (float): Resistenza dinamica Zener (resistenza nella regione di rottura).
        """
        super().__init__(name, anode_node, cathode_node)
        self.Is = Is
        self.N = N
        self.Vt = Vt
        self.Vz = Vz
        self.Iz = Iz
        self.Rz = Rz

    def calculate_current(self, Vd: float) -> float:
        """
        Calcola la corrente attraverso il diodo Zener data la tensione ai suoi capi (Vd = V_anode - V_cathode).
        Modella la conduzione in avanti e la rottura Zener in inversa.
        """
        current = 0.0

        # Conduzione in avanti (Vd > 0) - come un diodo normale
        if Vd >= 0:
            if Vd / (self.N * self.Vt) > 700:
                current = self.Is * (np.exp(700) - 1)
            else:
                current = self.Is * (np.exp(Vd / (self.N * self.Vt)) - 1)
        # Rottura Zener in inversa (Vd < 0 e |Vd| > Vz)
        else:
            # Modello semplificato per la rottura Zener
            # La corrente aumenta rapidamente una volta superata Vz
            if Vd < -self.Vz:
                # Modello lineare nella regione di rottura
                current = (Vd + self.Vz) / self.Rz - self.Iz # La corrente è negativa
            else:
                # Corrente di saturazione inversa (molto piccola) prima della rottura
                current = -self.Is # Corrente inversa negativa

        return current

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

