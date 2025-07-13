import numpy as np

class Resistor:
    """
    Modello fisico per un resistore.
    """
    def __init__(self, resistance=1e3):
        """
        Inizializza il resistore.

        Args:
            resistance (float): Il valore di resistenza in Ohm. Default 1 kOhm.
        """
        if resistance <= 0:
            raise ValueError("La resistenza deve essere un valore positivo.")
        self.R = float(resistance)
        print(f"Resistore creato con R = {self.R} Ohm.")

    def get_resistance(self):
        """
        Restituisce il valore di resistenza del componente.
        """
        return self.R

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il resistore data la differenza di potenziale.
        (Legge di Ohm: I = V/R)

        Args:
            voltage_difference (float or np.ndarray): Differenza di potenziale (V) ai capi del resistore.

        Returns:
            float or np.ndarray: Corrente (A) che scorre attraverso il resistore.
        """
        return voltage_difference / self.R

    def __str__(self):
        return f"Resistore(R={self.R} Ohm)"

# Esempio di utilizzo (puoi metterlo in main.py o in un test)
if __name__ == "__main__":
    r1 = Resistor(1000)
    v_diff = 5  # Volt
    current = r1.calculate_current(v_diff)
    print(f"Con {v_diff}V, la corrente è {current}A")

    r2 = Resistor(resistance=2.2e3) # 2.2 kOhm
    print(r2)

    try:
        r_invalid = Resistor(0)
    except ValueError as e:
        print(f"Errore previsto: {e}")
