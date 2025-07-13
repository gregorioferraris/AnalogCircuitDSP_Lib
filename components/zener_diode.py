import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM # Importa costanti fisiche

class ZenerDiode:
    """
    Modello fisico per un diodo Zener.
    Combina l'equazione di Shockley per la polarizzazione diretta
    e un modello per la regione di breakdown Zener in polarizzazione inversa.
    """
    def __init__(self, Is=1e-12, n=1.0, Vz=5.6, Iz_test=1e-3, Rz=5.0, temperature_k=T_ROOM):
        """
        Inizializza il diodo Zener.

        Args:
            Is (float): Corrente di saturazione inversa (A) per la polarizzazione diretta.
            n (float): Fattore di idealità per la polarizzazione diretta.
            Vz (float): Tensione Zener nominale (V).
            Iz_test (float): Corrente di test (A) alla quale è specificato Vz.
            Rz (float): Resistenza dinamica nella regione Zener (Ohm).
            temperature_k (float): Temperatura in Kelvin.
        """
        if Is <= 0 or n <= 0 or Vz <= 0 or Iz_test <= 0 or Rz < 0 or temperature_k <= 0:
            raise ValueError("I parametri Is, n, Vz, Iz_test, Rz, Temperature devono essere positivi o non negativi.")

        self.Is = float(Is)
        self.n = float(n)
        self.Vz = float(Vz)
        self.Iz_test = float(Iz_test)
        self.Rz = float(Rz)
        self.Vt = (K * temperature_k) / Q # Tensione termica

        # Calcola la tensione di knee (Vzk) per la transizione alla regione Zener
        # Usiamo un punto per definire l'inizio della curva Zener
        # Questo è un modello semplificato; modelli più precisi usano un esponenziale per la regione Zener.
        # Qui usiamo un'approssimazione lineare per la parte Zener per semplicità.
        self.Vzk = self.Vz - (self.Iz_test * self.Rz) # Tensione all'inizio della pendenza Zener

        print(f"Zener Diode creato con Vz={self.Vz}V, Rz={self.Rz} Ohm.")

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il diodo Zener data la differenza di potenziale.

        Args:
            voltage_difference (float or np.ndarray): Tensione (Vd) ai capi del diodo.

        Returns:
            float or np.ndarray: Corrente (Id) che scorre attraverso il diodo.
        """
        Id = np.zeros_like(voltage_difference, dtype=float) if isinstance(voltage_difference, np.ndarray) else 0.0

        # Regione di polarizzazione diretta (forward bias)
        forward_condition = voltage_difference > 0
        if np.any(forward_condition):
            Vd_forward = voltage_difference[forward_condition]
            exponent_arg = Vd_forward / (self.n * self.Vt)
            exponent_arg = np.clip(exponent_arg, None, 700)
            Id[forward_condition] = self.Is * (np.exp(exponent_arg) - 1)

        # Regione di polarizzazione inversa (reverse bias)
        reverse_condition = voltage_difference <= 0
        if np.any(reverse_condition):
            Vd_reverse = voltage_difference[reverse_condition]

            # Corrente inversa di saturazione (vicina a -Is)
            # Fino al punto in cui si raggiunge la tensione Zener
            zener_breakdown_condition = Vd_reverse <= -self.Vzk
            if np.any(zener_breakdown_condition):
                # Modello lineare nella regione Zener
                Id[reverse_condition & zener_breakdown_condition] = (Vd_reverse[reverse_condition & zener_breakdown_condition] + self.Vz) / self.Rz + self.Iz_test
            else:
                # Corrente di fuga inversa (quasi -Is)
                Id[reverse_condition & ~zener_breakdown_condition] = -self.Is

        return Id

    def __str__(self):
        return f"ZenerDiode(Vz={self.Vz}V, Iz_test={self.Iz_test}A, Rz={self.Rz}Ohm)"

# Esempio di utilizzo e plot delle curve I-V
if __name__ == "__main__":
    zener1 = ZenerDiode(Vz=5.6, Iz_test=1e-3, Rz=10.0) # Zener 5.6V
    print(zener1)

    import matplotlib.pyplot as plt

    voltages = np.linspace(-10, 1, 200) # Da -10V a 1V per mostrare Zener e forward
    currents = np.array([zener1.calculate_current(v) for v in voltages])

    plt.figure(figsize=(10, 6))
    plt.plot(voltages, currents * 1e3) # Corrente in mA
    plt.title(f'Caratteristica I-V del Diodo Zener ({zener1.Vz}V)')
    plt.xlabel('Tensione (V)')
    plt.ylabel('Corrente (mA)')
    plt.grid(True)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.annotate(f'Vz = {-zener1.Vz}V', xy=(-zener1.Vz, 0), xytext=(-zener1.Vz - 2, 5),
                 arrowprops=dict(facecolor='black', shrink=0.05),
                 horizontalalignment='right')
    plt.show()
