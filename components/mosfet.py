import numpy as np

class MOSFET:
    """
    Modello fisico semplificato per un MOSFET a canale N ad enhancement mode.
    Copre le regioni di Cutoff, Triodo (Ohmic) e Saturazione.
    """
    def __init__(self, Vt=1.0, Kn=0.5e-3, lambda_val=0.0):
        """
        Inizializza il MOSFET.

        Args:
            Vt (float): Tensione di soglia (Threshold Voltage) in Volt. Default 1.0V.
            Kn (float): Transconduttanza del dispositivo (k' * W/L) in A/V^2. Default 0.5mA/V^2.
            lambda_val (float): Parametro di modulazione della lunghezza del canale.
                                  (Per modellare la pendenza delle curve in saturazione). Default 0.0 (ideale).
        """
        if Vt < 0 and Kn > 0: # Enhancement mode (Vt > 0 for N-channel)
            raise ValueError("La tensione di soglia (Vt) deve essere positiva per MOSFET N-channel ad enhancement mode.")
        if Kn <= 0:
            raise ValueError("Kn deve essere un valore positivo.")
        if not (0 <= lambda_val < 1):
            raise ValueError("Lambda deve essere tra 0 e 1 (escluso 1).")

        self.Vt = float(Vt)
        self.Kn = float(Kn)
        self.lambda_val = float(lambda_val)

        print(f"MOSFET (N-Enh) creato con Vt={self.Vt}V, Kn={self.Kn}A/V^2, lambda={self.lambda_val}.")

    def calculate_drain_current(self, Vgs, Vds):
        """
        Calcola la corrente di Drain (Id) del MOSFET.

        Args:
            Vgs (float or np.ndarray): Tensione Gate-Source (V).
            Vds (float or np.ndarray): Tensione Drain-Source (V).

        Returns:
            float or np.ndarray: Corrente di Drain (A).
        """
        Id = np.zeros_like(Vgs, dtype=float) if isinstance(Vgs, np.ndarray) else 0.0

        # Overdrive voltage (tensione oltre la soglia)
        Vov = Vgs - self.Vt

        # Regione di Cutoff (Id = 0)
        # Il transistor è spento se Vgs <= Vt
        cutoff_condition = Vov <= 0
        Id[cutoff_condition] = 0.0

        # Regione di Saturazione
        # Il transistor è acceso e Vds >= Vov (o Vds >= Vgs - Vt)
        saturation_condition = ~cutoff_condition & (Vds >= Vov)
        if np.any(saturation_condition):
            Vov_sat = Vov[saturation_condition]
            Vds_sat = Vds[saturation_condition]
            Id[saturation_condition] = 0.5 * self.Kn * (Vov_sat**2) * (1 + self.lambda_val * Vds_sat)

        # Regione di Triodo (Ohmic)
        # Il transistor è acceso e Vds < Vov
        triode_condition = ~cutoff_condition & (Vds < Vov)
        if np.any(triode_condition):
            Vov_tri = Vov[triode_condition]
            Vds_tri = Vds[triode_condition]
            Id[triode_condition] = self.Kn * ((Vov_tri * Vds_tri) - (0.5 * Vds_tri**2)) * (1 + self.lambda_val * Vds_tri)

        # Assicurati che la corrente non sia mai negativa
        return np.maximum(0, Id)

    def __str__(self):
        return f"MOSFET(Vt={self.Vt}V, Kn={self.Kn:.2e}A/V^2, lambda={self.lambda_val})"

# Esempio di utilizzo e plot delle curve I-V
if __name__ == "__main__":
    mosfet1 = MOSFET(Vt=1.5, Kn=1.0e-3, lambda_val=0.02) # MOSFET N-Enh tipico
    print(mosfet1)

    import matplotlib.pyplot as plt

    # Plot delle curve Id vs Vds per diversi Vgs
    vds_values = np.linspace(0, 10, 100) # Vds da 0 a 10V
    vgs_steps = np.linspace(mosfet1.Vt, mosfet1.Vt + 5, 5) # Diversi valori di Vgs sopra la soglia

    plt.figure(figsize=(10, 6))
    for vgs in vgs_steps:
        id_values = mosfet1.calculate_drain_current(vgs, vds_values)
        plt.plot(vds_values, id_values * 1e3, label=f'Vgs={vgs:.1f}V') # Corrente in mA

    plt.title('Caratteristiche di Uscita del MOSFET (Id vs Vds)')
    plt.xlabel('Vds (V)')
    plt.ylabel('Id (mA)')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot delle curve di trasferimento (Id vs Vgs per Vds fissa)
    vgs_values_transfer = np.linspace(0, mosfet1.Vt + 3, 100)
    vds_fixed = 5.0 # Vds fisso in regione di saturazione
    id_transfer = mosfet1.calculate_drain_current(vgs_values_transfer, vds_fixed)

    plt.figure(figsize=(8, 5))
    plt.plot(vgs_values_transfer, id_transfer * 1e3)
    plt.title('Caratteristiche di Trasferimento del MOSFET (Id vs Vgs)')
    plt.xlabel('Vgs (V)')
    plt.ylabel('Id (mA)')
    plt.grid(True)
    plt.show()
