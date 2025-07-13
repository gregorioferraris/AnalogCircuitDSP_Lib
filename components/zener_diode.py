# components/zener_diode.py

import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM

class ZenerDiode:
    def __init__(self, Vz=5.6, Iz_test=1e-3, Rz=10.0, saturation_current=1e-14, emission_coefficient=1.0):
        """
        Simula un diodo Zener.
        Args:
            Vz (float): Tensione Zener nominale (Volt).
            Iz_test (float): Corrente di test a cui Vz è specificata (Ampere).
            Rz (float): Resistenza dinamica nella regione Zener (Ohm).
            saturation_current (float): Corrente di saturazione inversa (per la regione forward).
            emission_coefficient (float): Coefficiente di emissione (per la regione forward).
        """
        self.Vz = float(Vz) # Zener voltage (magnitudine)
        self.Iz_test = float(Iz_test) # Test current for Zener voltage
        self.Rz = float(Rz) # Zener dynamic resistance
        self.Is = float(saturation_current) # Forward saturation current
        self.n = float(emission_coefficient) # Forward emission coefficient
        self.Vt = VT_ROOM # Thermal voltage

        # Calcola un parametro Vbd per la regione di breakdown Zener
        # Vz = Vbd + Iz_test * Rz
        # Vbd = Vz - Iz_test * Rz (approssimazione)
        # Usiamo un modello basato sul diodo per la regione Zener inversa.
        # Un modo comune è modellare la regione inversa come:
        # I_R = -Iz * exp(-(V_z + Vd) / (n_z * Vt))
        # dove V_z è la tensione zener e Vd è la tensione applicata (negativa).
        # Per semplicità, useremo un modello di diodo inverso con una resistenza.

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il diodo Zener.
        Considera la regione di forward e la regione di breakdown Zener.
        """
        if voltage_difference >= 0:
            # Regione di polarizzazione diretta (come un diodo normale)
            exponent_arg = voltage_difference / (self.n * self.Vt)
            if isinstance(exponent_arg, np.ndarray):
                exponent_arg = np.clip(exponent_arg, None, 700)
            else:
                exponent_arg = min(exponent_arg, 700)
            return self.Is * (np.exp(exponent_arg) - 1.0)
        else:
            # Regione di polarizzazione inversa (Zener)
            # La tensione Zener è Vz, ma la tensione applicata è negativa (-V_app).
            # Se |V_app| > Vz, la corrente inversa aumenta rapidamente.
            # Modello semplificato: una volta superato Vz, la corrente è guidata dalla resistenza Rz.
            if voltage_difference < -self.Vz: # Entrato nella regione Zener
                # Corrente = (Tensione applicata - Tensione Zener) / Resistenza Zener
                # La tensione applicata è negativa, quindi V_diff - (-Vz)
                return (voltage_difference + self.Vz) / self.Rz # Corrente negativa
            else:
                # Corrente di fuga molto piccola in polarizzazione inversa prima del breakdown
                return -self.Is # -Is è la corrente di fuga inversa

    def calculate_conductance(self, voltage_difference):
        """
        Calcola la conduttanza dinamica (differenziale) del diodo Zener.
        """
        if voltage_difference >= 0:
            # Conduttanza in polarizzazione diretta (come diodo normale)
            exponent_arg = voltage_difference / (self.n * self.Vt)
            if isinstance(exponent_arg, np.ndarray):
                exponent_arg = np.clip(exponent_arg, None, 700)
            else:
                exponent_arg = min(exponent_arg, 700)
            return (self.Is / (self.n * self.Vt)) * np.exp(exponent_arg)
        else:
            # Conduttanza in polarizzazione inversa (Zener)
            if voltage_difference < -self.Vz:
                return 1.0 / self.Rz # Conduttanza dinamica nella regione Zener
            else:
                return 1e-9 # Conduttanza molto bassa (quasi zero) in regione di fuga inversa

    def __str__(self):
        return f"ZenerDiode(Vz={self.Vz:.1f}V, Rz={self.Rz:.1f} Ohm)"
