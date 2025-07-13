import numpy as np

class JFET:
    """
    Modello fisico per un JFET (Junction Field-Effect Transistor).
    Modello semplice che copre le regioni di Cutoff, Triodo (Ohmic) e Saturazione.
    """
    def __init__(self, Idss=0.01, Vp=-2.0, lambda_val=0.0):
        """
        Inizializza il JFET.

        Args:
            Idss (float): Corrente di drain a Vgs=0 e saturazione (A). Default 10mA.
            Vp (float): Tensione di pinch-off (Vgs a cui Id diventa quasi zero in saturazione) (V).
                        Solitamente negativo per N-channel JFET. Default -2.0V.
            lambda_val (float): Parametro di modulazione della lunghezza del canale.
                                  (Per modellare la pendenza delle curve in saturazione). Default 0.0 (ideale).
        """
        if Idss <= 0:
            raise ValueError("IDSS deve essere un valore positivo.")
        if Vp >= 0 and Idss > 0: # Vp < 0 for N-channel JFETs, Vp > 0 for P-channel JFETs (Idss > 0)
            raise ValueError("VP deve essere negativo per JFET N-channel con IDSS positivo (o viceversa).")
        if not (0 <= lambda_val < 1): # Common range for lambda
            raise ValueError("Lambda deve essere tra 0 e 1 (escluso 1).")

        self.Idss = float(Idss)
        self.Vp = float(Vp)
        self.lambda_val = float(lambda_val)

        print(f"JFET creato con Idss={self.Idss}A, Vp={self.Vp}V, lambda={self.lambda_val}.")

    def calculate_drain_current(self, Vgs, Vds):
        """
        Calcola la corrente di Drain (Id) del JFET.

        Args:
            Vgs (float or np.ndarray): Tensione Gate-Source (V).
            Vds (float or np.ndarray): Tensione Drain-Source (V).

        Returns:
            float or np.ndarray: Corrente di Drain (A).
        """
        Id = np.zeros_like(Vgs, dtype=float) if isinstance(Vgs, np.ndarray) else 0.0

        # Clipping Vgs per evitare problemi numerici e rispettare i limiti fisici
        # Il gate di un JFET non dovrebbe mai andare positivo rispetto al source oltre un certo limite
        # per non entrare in conduzione la giunzione gate-source.
        # Per semplicità, consideriamo Vgs <= 0 per la conduzione normale del canale.
        Vgs_effective = np.minimum(Vgs, 0.0) # Assume JFET N-channel

        # Regione di Cutoff (Id = 0)
        # Il canale è chiuso se Vgs <= Vp
        cutoff_condition = Vgs_effective <= self.Vp
        Id[cutoff_condition] = 0.0

        # Regione di Saturazione (Id = Idss * (1 - Vgs/Vp)^2 * (1 + lambda * Vds))
        # Vds deve essere sufficientemente grande (Vds >= Vgs - Vp)
        saturation_condition = ~cutoff_condition & (Vds >= (Vgs_effective - self.Vp))
        if np.any(saturation_condition):
            Vgs_norm = (1 - Vgs_effective[saturation_condition] / self.Vp)
            Id[saturation_condition] = self.Idss * (Vgs_norm**2) * (1 + self.lambda_val * Vds[saturation_condition])

        # Regione di Triodo (Ohmic)
        # Vds è più piccola (Vds < Vgs - Vp) e il canale è aperto
        triode_condition = ~cutoff_condition & (Vds < (Vgs_effective - self.Vp))
        if np.any(triode_condition):
            Vds_eff = Vds[triode_condition]
            Vgs_eff = Vgs_effective[triode_condition]
            Id[triode_condition] = self.Idss * (2 * (1 - Vgs_eff / self.Vp) * (Vds_eff / self.Vp) - (Vds_eff / self.Vp)**2) * (1 + self.lambda_val * Vds_eff)

        # Assicurati che la corrente non sia mai negativa (non può scorrere al contrario nel canale)
        return np.maximum(0, Id)

    def __str__(self):
        return f"JFET(Idss={self.Idss}A, Vp={self.Vp}V, lambda={self.lambda_val})"

# Esempio di utilizzo e plot delle curve I-V
if __name__ == "__main__":
    jfet1 = JFET(Idss=0.01, Vp=-2.0, lambda_val=0.05) # JFET tipico N-channel
    print(jfet1)

    import matplotlib.pyplot as plt

    # Plot delle curve Id vs Vds per diversi Vgs
    vds_values = np.linspace(0, 5, 100) # Vds da 0 a 5V
    vgs_steps = np.linspace(jfet1.Vp, 0, 5) # Diversi valori di Vgs

    plt.figure(figsize=(10, 6))
    for vgs in vgs_steps:
        id_values = jfet1.calculate_drain_current(vgs, vds_values)
        plt.plot(vds_values, id_values * 1e3, label=f'Vgs={vgs:.1f}V') # Corrente in mA

    plt.title('Caratteristiche di Uscita del JFET (Id vs Vds)')
    plt.xlabel('Vds (V)')
    plt.ylabel('Id (mA)')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot delle curve di trasferimento (Id vs Vgs per Vds fissa)
    vgs_values_transfer = np.linspace(jfet1.Vp - 0.5, 0.0, 100)
    vds_fixed = 5.0 # Vds fisso in regione di saturazione
    id_transfer = jfet1.calculate_drain_current(vgs_values_transfer, vds_fixed)

    plt.figure(figsize=(8, 5))
    plt.plot(vgs_values_transfer, id_transfer * 1e3)
    plt.title('Caratteristiche di Trasferimento del JFET (Id vs Vgs)')
    plt.xlabel('Vgs (V)')
    plt.ylabel('Id (mA)')
    plt.grid(True)
    plt.show()
