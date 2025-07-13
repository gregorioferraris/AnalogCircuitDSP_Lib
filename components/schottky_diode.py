# components/schottky_diode.py

import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM

class SchottkyDiode:
    def __init__(self, saturation_current=1e-10, emission_coefficient=1.05):
        """
        Simula un diodo Schottky.
        Ha una caduta di tensione diretta inferiore e una riaccensione più veloce.
        Args:
            saturation_current (float): Corrente di saturazione inversa (più alta dei diodi PN).
            emission_coefficient (float): Coefficiente di emissione (vicino a 1).
        """
        self.Is = float(saturation_current)
        self.n = float(emission_coefficient)
        self.Vt = VT_ROOM

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il diodo Schottky usando il modello di Shockley.
        """
        # Limita l'esponente per prevenire overflow con tensioni elevate
        exponent_arg = voltage_difference / (self.n * self.Vt)
        if isinstance(exponent_arg, np.ndarray):
            exponent_arg = np.clip(exponent_arg, None, 700) # Limite ragionevole per np.exp
        else:
            exponent_arg = min(exponent_arg, 700) # Limite per singolo valore

        return self.Is * (np.exp(exponent_arg) - 1.0)

    def calculate_conductance(self, voltage_difference):
        """
        Calcola la conduttanza dinamica (differenziale) del diodo Schottky.
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
        return f"SchottkyDiode(Is={self.Is:.1e}, n={self.n:.2f})"
