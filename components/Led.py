# components/led.py

import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM

class LED:
    def __init__(self, saturation_current=1e-20, emission_coefficient=3.4, breakdown_voltage=-5.0, luminous_efficiency=100.0):
        # Parametri del diodo (diodo con una caduta di tensione maggiore)
        self.Is = float(saturation_current) # Corrente di saturazione inversa (molto bassa per LED)
        self.n = float(emission_coefficient) # Coefficiente di emissione più alto per LED
        self.Vt = VT_ROOM # Tensione termica
        self.Vf_typical = 1.8 # Tensione di forward tipica per l'accensione
        self.Rs = 1.0 # Resistenza in serie interna per modellare la pendenza dopo Vf

        # Parametri specifici del LED
        self.breakdown_voltage = float(breakdown_voltage) # Tensione di breakdown inverso
        self.luminous_efficiency = float(luminous_efficiency) # Efficienza luminosa (un fattore arbitrario)

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il LED.
        Include un modello semplificato per la caduta diretta e il breakdown inverso.
        """
        if voltage_difference > 0:
            # Modello del diodo con resistenza in serie per la regione di conduzione
            # Risolvi per Vd: Vd = V - I*Rs
            # Id = Is * (exp(Vd / (n*Vt)) - 1)
            # Questo è un'equazione implicita. Possiamo usare Newton-Raphson internamente
            # o una semplificazione. Per semplicità e compatibilità con MNA, usiamo una forma
            # che assomiglia a un diodo in serie con una resistenza.
            # Per l'MNA, il contributo sarà già incluso nella matrice Jacobiana
            # per la linearizzazione del diodo.
            # La resistenza in serie può essere modellata come un componente separato.
            # Qui si simula solo la corrente di un diodo senza Rs esplicita:
            exponent_arg = voltage_difference / (self.n * self.Vt)
            if isinstance(exponent_arg, np.ndarray):
                exponent_arg = np.clip(exponent_arg, None, 700)
            else:
                exponent_arg = min(exponent_arg, 700)
            return self.Is * (np.exp(exponent_arg) - 1.0)
        elif voltage_difference < self.breakdown_voltage:
            # Regione di breakdown inverso (corrente negativa crescente)
            # Modello semplificato: aumento rapido della corrente oltre la tensione di breakdown
            # Non è un modello Zener completo, ma una risposta inversa.
            return (voltage_difference - self.breakdown_voltage) / 10.0 # Rinv = 10 Ohm
        else:
            return 0.0 # Quasi zero in polarizzazione inversa e vicino a zero in forward prima di Vf

    def calculate_conductance(self, voltage_difference):
        """
        Calcola la conduttanza dinamica del LED.
        La derivata di I = Is * (exp(Vd/(n*Vt)) - 1) è g = Is / (n*Vt) * exp(Vd/(n*Vt))
        Aggiungiamo la derivata per la regione di breakdown.
        """
        if voltage_difference > 0:
            exponent_arg = voltage_difference / (self.n * self.Vt)
            if isinstance(exponent_arg, np.ndarray):
                exponent_arg = np.clip(exponent_arg, None, 700)
            else:
                exponent_arg = min(exponent_arg, 700)
            return (self.Is / (self.n * self.Vt)) * np.exp(exponent_arg)
        elif voltage_difference < self.breakdown_voltage:
            return 1.0 / 10.0 # Derivata di (V - V_bd) / 10 è 1/10
        else:
            return 1e-9 # Conduttanza molto bassa (quasi zero) quando è spento

    def get_luminous_output(self, current):
        """
        Calcola l'output luminoso relativo basato sulla corrente in forward.
        """
        # La luminosità è approssimativamente proporzionale alla corrente.
        # Solo per correnti positive.
        if isinstance(current, np.ndarray):
            luminous = np.maximum(0, current) * self.luminous_efficiency
        else:
            luminous = max(0, current) * self.luminous_efficiency
        return luminous

    def __str__(self):
        return f"LED(Is={self.Is:.1e}, n={self.n:.1f}, Vf_breakdown={self.breakdown_voltage:.1f}V)"
