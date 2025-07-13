import numpy as np

class RectifierTube:
    """
    Modello fisico semplificato per un Diodo Raddrizzatore a Valvola (Rectifier Tube).
    Modella principalmente la caduta di tensione in avanti e la resistenza interna.
    Non include modellazione complessa di tempo di riscaldamento o limitazioni di corrente dinamiche.
    """
    def __init__(self, forward_voltage_drop=15.0, series_resistance=50.0, saturation_current=1e-6):
        """
        Inizializza il diodo raddrizzatore a valvola.

        Args:
            forward_voltage_drop (float): Caduta di tensione tipica in avanti (V).
                                           Questo è un valore approssimativo per il "knee".
            series_resistance (float): Resistenza interna equivalente della valvola in conduzione (Ohm).
            saturation_current (float): Corrente di fuga inversa (A).
        """
        if forward_voltage_drop < 0 or series_resistance <= 0 or saturation_current <= 0:
            raise ValueError("I parametri devono essere positivi.")

        self.Vf = float(forward_voltage_drop)
        self.Rs = float(series_resistance)
        self.Is = float(saturation_current) # Per modellare la corrente inversa

        print(f"Rectifier Tube creato con Vf={self.Vf}V, Rs={self.Rs}Ohm.")

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso il diodo raddrizzatore.
        Modello piecewise linear o esponenziale per la regione di conduzione.
        Qui usiamo un modello semplice: spento, acceso con Vf e Rs.

        Args:
            voltage_difference (float or np.ndarray): Tensione (Vd) ai capi del diodo.

        Returns:
            float or np.ndarray: Corrente (Id) che scorre attraverso il diodo.
        """
        Id = np.zeros_like(voltage_difference, dtype=float) if isinstance(voltage_difference, np.ndarray) else 0.0

        # Regione di polarizzazione diretta
        forward_condition = voltage_difference >= self.Vf
        if np.any(forward_condition):
            # Corrente = (Tensione_oltre_Vf) / Rs
            Id[forward_condition] = (voltage_difference[forward_condition] - self.Vf) / self.Rs

        # Regione di blocco inverso o sotto soglia
        reverse_condition = voltage_difference < self.Vf
        if np.any(reverse_condition):
            # Per tensioni negative, c'è solo una piccola corrente di fuga inversa
            Id[reverse_condition] = -self.Is # Semplificato, un diodo ideale non conducerebbe affatto

        return Id

    def __str__(self):
        return f"RectifierTube(Vf={self.Vf}V, Rs={self.Rs}Ohm)"

# Esempio di utilizzo e plot delle curve I-V
if __name__ == "__main__":
    rectifier1 = RectifierTube(forward_voltage_drop=15.0, series_resistance=100.0) # Tipico 5AR4
    print(rectifier1)

    import matplotlib.pyplot as plt

    voltages = np.linspace(-50, 50, 100) # Da -50V a 50V
    currents = np.array([rectifier1.calculate_current(v) for v in voltages])

    plt.figure(figsize=(10, 6))
    plt.plot(voltages, currents * 1e3) # Corrente in mA
    plt.title('Caratteristica I-V del Diodo Raddrizzatore a Valvola')
    plt.xlabel('Tensione (V)')
    plt.ylabel('Corrente (mA)')
    plt.grid(True)
    plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    plt.show()
