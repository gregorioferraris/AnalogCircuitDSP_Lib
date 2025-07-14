# components/pentode.py
import numpy as np
from components.component import Component

class Pentode(Component):
    def __init__(self, name: str, plate_node: str, grid_node: str, screen_grid_node: str, suppressor_grid_node: str, cathode_node: str,
                 mu: float = 100.0, Kp: float = 1.0, Ex: float = 1.0, Kg1: float = 1.0, Kg2: float = 1.0, Vg_offset: float = 0.0):
        """
        Inizializza un Pentodo usando un modello tipo Koren esteso (o simile).
        Args:
            name (str): Nome univoco dell'istanza.
            plate_node (str): Nodo della placca (anodo).
            grid_node (str): Nodo della griglia di controllo (g1).
            screen_grid_node (str): Nodo della griglia schermo (g2).
            suppressor_grid_node (str): Nodo della griglia soppressore (g3).
            cathode_node (str): Nodo del catodo.
            mu (float): Fattore di amplificazione (per g1).
            Kp (float): Parametro di permeabilità.
            Ex (float): Esponente.
            Kg1 (float): Parametro di linearità della griglia di controllo.
            Kg2 (float): Parametro di linearità della griglia schermo.
            Vg_offset (float): Offset di tensione per la griglia di controllo.
        """
        super().__init__(name, plate_node, grid_node, screen_grid_node, suppressor_grid_node, cathode_node)
        self.pin_names = ('plate', 'grid', 'screen_grid', 'suppressor_grid', 'cathode') # Per mappatura nodi con nome
        self.mu = mu
        self.Kp = Kp
        self.Ex = Ex
        self.Kg1 = Kg1
        self.Kg2 = Kg2
        self.Vg_offset = Vg_offset

    def calculate_plate_current(self, Vgk: float, Vpk: float, Vg2k: float, Vg3k: float) -> float:
        """
        Calcola la corrente di placca (Ip) per un pentodo.
        Questo è un modello molto semplificato; i modelli pentodo reali sono complessi.
        Args:
            Vgk (float): Tensione Griglia-Catodo.
            Vpk (float): Tensione Placca-Catodo.
            Vg2k (float): Tensione Griglia Schermo-Catodo.
            Vg3k (float): Tensione Griglia Soppressore-Catodo.
        Returns:
            float: Corrente di placca in Ampere.
        """
        # Per semplicità, questo modello ignora Vg3k e usa un approccio simile al triodo,
        # ma con l'influenza della griglia schermo.
        # Un modello pentodo più accurato richiederebbe equazioni più complesse
        # che modellano la partizione di corrente tra placca e griglia schermo.

        if Vpk <= 0:
            return 0.0 # Nessuna corrente se la placca non è polarizzata positivamente

        # Calcolo della tensione di griglia equivalente che controlla la corrente totale
        # Questo è un punto di semplificazione, i modelli reali sono più complessi
        effective_grid_voltage = (Vgk + self.Vg_offset) + (Vg2k / self.mu)

        if effective_grid_voltage <= 0:
            return 0.0 # Cutoff

        # La corrente di placca è influenzata dalla tensione di placca e di griglia schermo
        # Questo è un modello comportamentale, non un modello fisico completo di pentodo.
        # Per un modello più rigoroso, si dovrebbe considerare la partizione di corrente
        # tra placca e griglia schermo e la saturazione di placca.
        
        # Un approccio comune è usare una funzione di saturazione per la placca (es. arctan, tanh)
        # per modellare il "ginocchio" della curva di placca.
        # Per ora, usiamo una dipendenza lineare dalla tensione di placca per semplicità,
        # ma questo è un punto da migliorare per fedeltà.
        
        # La corrente totale (placca + griglia schermo) dipende principalmente da Vgk e Vg2k
        total_current = self.Kp * (effective_grid_voltage)**self.Ex

        # Partizione della corrente (molto semplificata: la corrente di placca è una frazione della totale)
        # In un pentodo, la corrente di placca è quasi indipendente da Vpk una volta in saturazione.
        # Qui, usiamo un fattore che dipende da Vpk per simulare una leggera dipendenza.
        # Questo è un placeholder e andrebbe sostituito con un modello pentodo più specifico (es. Koren Pentode, o modelli basati su curve).
        plate_current_factor = np.tanh(Vpk / 50.0) # Esempio di funzione di saturazione per Vpk

        Ip = total_current * plate_current_factor
        
        # Aggiungi un limite per evitare valori irrealistici
        if Ip > 1.0:
            Ip = 1.0

        return Ip

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

