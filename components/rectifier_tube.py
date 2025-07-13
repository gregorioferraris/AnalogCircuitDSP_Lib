# components/rectifier_tube.py

import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM # Importa costanti fisiche

class RectifierTube:
    def __init__(self, forward_voltage_drop=15.0, series_resistance=50.0, saturation_current=1e-6, emission_coefficient=2.0):
        """
        Simula un diodo raddrizzatore a vuoto (es. 5Y3, 5AR4).
        Questi diodi hanno una caduta di tensione diretta significativa e una certa resistenza interna.
        Args:
            forward_voltage_drop (float): Caduta di tensione tipica in forward (Volt).
            series_resistance (float): Resistenza interna in serie in conduzione (Ohm).
            saturation_current (float): Corrente di saturazione inversa (più alta dei semiconduttori).
            emission_coefficient (float): Coefficiente di emissione (più alto dei semiconduttori).
        """
        self.Vf = float(forward_voltage_drop) # Tipica caduta di tensione in conduzione
        self.Rs = float(series_resistance) # Resistenza in serie interna
        self.Is = float(saturation_current) # Corrente di saturazione inversa
        self.n = float(emission_coefficient) # Coefficiente di emissione (alto per modellare Vf)
        self.Vt = VT_ROOM # Tensione termica

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il diodo raddrizzatore.
        Modello semplificato che combina una caduta di tensione fissa e una resistenza in serie,
        o un modello di diodo più generale.
        """
        if voltage_difference > self.Vf: # Se la tensione supera la caduta
            # Regione di conduzione: (V - Vf) / Rs
            # Questo è un modello lineare dopo la caduta di tensione.
            return (voltage_difference - self.Vf) / self.Rs
        elif voltage_difference > 0:
            # Tra 0 e Vf, una curva graduale (approssimazione di diodo)
            # Usa il modello diodo standard per questa transizione.
            exponent_arg = voltage_difference / (self.n * self.Vt)
            if isinstance(exponent_arg, np.ndarray):
                exponent_arg = np.clip(exponent_arg, None, 700)
            else:
                exponent_arg = min(exponent_arg, 700)
            return self.Is * (np.exp(exponent_arg) - 1.0)
        else:
            # Regione di blocco inverso o quasi zero
            return -self.Is # Corrente di fuga inversa (Is è positiva)

    def calculate_conductance(self, voltage_difference):
        """
        Calcola la conduttanza dinamica (differenziale) del diodo raddrizzatore.
        """
        if voltage_difference > self.Vf:
            # Conduttanza nella regione di conduzione: 1/Rs
            return 1.0 / self.Rs
        elif voltage_difference > 0:
            # Conduttanza nella regione di transizione (diodo)
            exponent_arg = voltage_difference / (self.n * self.Vt)
            if isinstance(exponent_arg, np.ndarray):
                exponent_arg = np.clip(exponent_arg, None, 700)
            else:
                exponent_arg = min(exponent_arg, 700)
            return (self.Is / (self.n * self.Vt)) * np.exp(exponent_arg)
        else:
            # Conduttanza molto bassa (quasi zero) in polarizzazione inversa
            return 1e-9 # Una conduttanza molto piccola per stabilità numerica

    def __str__(self):
        return f"RectifierTube(Vf={self.Vf:.1f}V, Rs={self.Rs:.1f} Ohm)"
