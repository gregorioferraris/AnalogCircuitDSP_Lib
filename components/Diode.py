import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM # Importa costanti fisiche

class Diode:
    """
    Modello fisico per un diodo standard (equazione di Shockley).
    """
    def __init__(self, saturation_current=1e-12, ideality_factor=1.0, temperature_k=T_ROOM):
        """
        Inizializza il diodo.

        Args:
            saturation_current (float): Corrente di saturazione inversa (Is) in Ampere. Default 1e-12A.
            ideality_factor (float): Fattore di idealità (n). Default 1.0.
            temperature_k (float): Temperatura in Kelvin. Default temperatura ambiente.
        """
        if saturation_current <= 0:
            raise ValueError("La corrente di saturazione deve essere un valore positivo.")
        if ideality_factor <= 0:
            raise ValueError("Il fattore di idealità deve essere un valore positivo.")
        if temperature_k <= 0:
            raise ValueError("La temperatura in Kelvin deve essere positiva.")

        self.Is = float(saturation_current)
        self.n = float(ideality_factor)
        self.Vt = (K * temperature_k) / Q # Tensione termica (kT/q)

        print(f"Diodo creato con Is = {self.Is}A, n = {self.n}, Vt = {self.Vt:.4f}V.")

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il diodo data la differenza di potenziale.
        (Equazione di Shockley: Id = Is * (exp(Vd / (n*Vt)) - 1))

        Args:
            voltage_difference (float or np.ndarray): Tensione (Vd) ai capi del diodo.

        Returns:
            float or np.ndarray: Corrente (Id) che scorre attraverso il diodo.
        """
        # Protezione per valori molto grandi di Vd per evitare overflow con np.exp
        # Questo è importante in contesti di solutori numerici.
        # Un valore ragionevole per Vd/(n*Vt) che causa overflow è tipicamente > 700
        exponent_arg = voltage_difference / (self.n * self.Vt)
        if isinstance(exponent_arg, np.ndarray):
            exponent_arg = np.clip(exponent_arg, None, 700) # Limita l'esponente per array
        else:
            exponent_arg = min(exponent_arg, 700) # Limita l'esponente per singolo valore

        return self.Is * (np.exp(exponent_arg) - 1)

    def __str__(self):
        return f"Diode(Is={self.Is:.2e}A, n={self.n}, Vt={self.Vt:.4f}V)"

# Esempio di utilizzo
if __name__ == "__main__":
    d1 = Diode() # Diodo standard con parametri default
    print(d1)

    voltages = np.linspace(-1, 1, 100) # Da -1V a 1V
    currents = [d1.calculate_current(v) for v in voltages]

    import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 5))
    plt.plot(voltages, currents)
    plt.title('Caratteristica I-V del Diodo')
    plt.xlabel('Tensione (V)')
    plt.ylabel('Corrente (A)')
    plt.grid(True)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.show()

    # Esempio di clipping
    signal = 0.5 * np.sin(np.linspace(0, 2*np.pi, 1000)) * 5 # Onda sinusoidale da -5V a 5V
    clipped_signal = [d1.calculate_current(s) for s in signal]
    
    plt.figure(figsize=(8, 5))
    plt.plot(signal, label='Input Signal (V)')
    plt.plot(clipped_signal, label='Output Current (A) - Clipped')
    plt.title('Esempio di Clipping del Diodo')
    plt.xlabel('Sample')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True)
    plt.show()
