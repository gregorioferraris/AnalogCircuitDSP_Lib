# components/diode.py

import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM # Importa costanti fisiche

class Diode:
    def __init__(self, saturation_current=1e-14, emission_coefficient=1.0):
        self.Is = float(saturation_current) # Corrente di saturazione inversa
        self.n = float(emission_coefficient) # Coefficiente di emissione (idealità)
        self.Vt = VT_ROOM # Tensione termica (kT/q)

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il diodo usando il modello di Shockley.
        Id = Is * (exp(Vd / (n*Vt)) - 1)
        """
        # Limita l'esponente per prevenire overflow con tensioni elevate
        exponent_arg = voltage_difference / (self.n * self.Vt)
        if isinstance(exponent_arg, np.ndarray):
            # Per gestire array di tensioni (es. per plotting)
            exponent_arg = np.clip(exponent_arg, None, 700) # Limite ragionevole per np.exp
        else:
            exponent_arg = min(exponent_arg, 700) # Limite per singolo valore

        return self.Is * (np.exp(exponent_arg) - 1.0)

    def calculate_conductance(self, voltage_difference):
        """
        Calcola la conduttanza dinamica (differenziale) del diodo.
        gd = d(Id)/d(Vd) = Is / (n*Vt) * exp(Vd / (n*Vt))
        """
        # Limita l'esponente per prevenire overflow
        exponent_arg = voltage_difference / (self.n * self.Vt)
        if isinstance(exponent_arg, np.ndarray):
            exponent_arg = np.clip(exponent_arg, None, 700)
        else:
            exponent_arg = min(exponent_arg, 700)

        return (self.Is / (self.n * self.Vt)) * np.exp(exponent_arg)

    def __str__(self):
        return f"Diode(Is={self.Is:.1e}, n={self.n:.1f})"
