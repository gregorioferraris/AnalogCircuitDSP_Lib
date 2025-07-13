import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM # Importa costanti fisiche

class SchottkyDiode:
    """
    Modello fisico per un diodo Schottky.
    Simile al diodo standard ma con una minore caduta di tensione e Is più alta.
    """
    def __init__(self, saturation_current=1e-6, ideality_factor=1.05, temperature_k=T_ROOM):
        """
        Inizializza il diodo Schottky.

        Args:
            saturation_current (float): Corrente di saturazione inversa (Is) in Ampere. Default 1e-6A (più alta).
            ideality_factor (float): Fattore di idealità (n). Default 1.05 (vicino a 1).
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

        print(f"Schottky Diode creato con Is = {self.Is}A, n = {self.n}, Vt = {self.Vt:.4f}V.")

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il diodo Schottky data la differenza di potenziale.
        (Equazione di Shockley: Id = Is * (exp(Vd / (n*Vt)) - 1))

        Args:
            voltage_difference (float or np.ndarray): Tensione (Vd) ai capi del diodo.

        Returns:
            float or np.ndarray: Corrente (Id) che scorre attraverso il diodo.
        """
        exponent_arg = voltage_difference / (self.n * self.Vt)
        if isinstance(exponent_arg, np.ndarray):
            exponent_arg = np.clip(exponent_arg, None, 700)
        else:
            exponent_arg = min(exponent_arg, 700)

        return self.Is * (np.exp(exponent_arg) - 1)

    def __str__(self):
        return f"SchottkyDiode(Is={self.Is:.2e}A, n={self.n}, Vt={self.Vt:.4f}V)"

# Esempio di utilizzo e plot delle curve I-V (confronto con diodo standard)
if __name__ == "__main__":
    schottky1 = SchottkyDiode()
    std_diode = Diode() # Diodo standard per confronto

    print(schottky1)

    voltages = np.linspace(-0.5, 1.0, 100) # Da -0.5V a 1.0V
    schottky_currents = np.array([schottky1.calculate_current(v) for v in voltages])
    std_diode_currents = np.array([std_diode.calculate_current(v) for v in voltages])

    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(voltages, schottky_currents * 1e3, label='Schottky Diode Current (mA)')
    plt.plot(voltages, std_diode_currents * 1e3, label='Standard Diode Current (mA)', linestyle='--')
    plt.title('Caratteristica I-V: Schottky vs Standard Diode')
    plt.xlabel('Tensione (V)')
    plt.ylabel('Corrente (mA)')
    plt.legend()
    plt.grid(True)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.show()
