# components/bjt.py
import numpy as np
from components.component import Component

class BJT(Component):
    def __init__(self, name: str, collector_node: str, base_node: str, emitter_node: str,
                 Is: float = 1e-14, Bf: float = 100.0, Br: float = 1.0, Vt: float = 0.0258, type: str = 'npn'):
        """
        Inizializza un BJT (Transistor a Giunzione Bipolare) usando un modello Ebers-Moll semplificato.
        Args:
            name (str): Nome univoco dell'istanza (es. "Q1").
            collector_node (str): Nome del nodo del collettore.
            base_node (str): Nome del nodo della base.
            emitter_node (str): Nome del nodo dell'emettitore.
            Is (float): Corrente di saturazione inversa.
            Bf (float): Guadagno di corrente in polarizzazione diretta (Beta Forward).
            Br (float): Guadagno di corrente in polarizzazione inversa (Beta Reverse).
            Vt (float): Tensione termica (kT/q).
            type (str): Tipo di BJT ('npn' o 'pnp').
        """
        super().__init__(name, collector_node, base_node, emitter_node)
        self.pin_names = ('collector', 'base', 'emitter') # Per mappatura nodi con nome
        self.Is = Is
        self.Bf = Bf
        self.Br = Br
        self.Vt = Vt
        self.type = type.lower()
        if self.type not in ['npn', 'pnp']:
            raise ValueError("BJT type must be 'npn' or 'pnp'.")

    def calculate_currents(self, Vbe: float, Vbc: float) -> tuple:
        """
        Calcola le correnti di collettore (Ic) e di base (Ib) per un BJT.
        Usa un modello Ebers-Moll semplificato.
        Args:
            Vbe (float): Tensione Base-Emettitore.
            Vbc (float): Tensione Base-Collettore.
        Returns:
            tuple: (Ic, Ib) - Corrente di collettore e corrente di base.
        """
        # Equazioni di Ebers-Moll semplificate
        # Corrente diodo emettitore-base
        Ife = self.Is * (np.exp(Vbe / self.Vt) - 1)
        # Corrente diodo collettore-base
        Irc = self.Is * (np.exp(Vbc / self.Vt) - 1)

        # Corrente di collettore
        Ic = (Ife / self.Bf) - Irc * (1 + 1 / self.Br)
        if self.type == 'pnp':
            Ic = -Ic # Inverti la direzione per PNP

        # Corrente di base
        Ib = (Ife / self.Bf) + (Irc / self.Br)
        if self.type == 'pnp':
            Ib = -Ib # Inverti la direzione per PNP

        # Corrente di emettitore (Ic + Ib)
        # Ie = Ife * (1 + 1 / self.Bf) - Irc * (1 + 1 / self.Br)

        return Ic, Ib

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i componenti non lineari, get_stamps restituisce matrici/vettori vuoti.
        Il loro contributo è gestito direttamente nel _system_equations del solutore.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

