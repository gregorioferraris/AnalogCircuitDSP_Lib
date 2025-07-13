import numpy as np

class Pentode:
    """
    Modello fisico per un Pentodo/Tetrodo (valvola termoionica a 5 o 4 elettrodi)
    basato su un'estensione del modello di Koren o un modello empirico comune.
    Per semplicità, iniziamo con un modello che usa Va, Vg1 e Vg2.
    """
    def __init__(self, mu=10.0, Kp=300.0, X=1.5, Kg1=5.0, Kg2=10.0, KvB=300.0):
        """
        Inizializza il Pentodo/Tetrodo con parametri estesi.

        Args:
            mu (float): Fattore di amplificazione (mu).
            Kp (float): Parametro di transconduttanza (mA/V^X).
            X (float): Esponente non lineare.
            Kg1 (float): Parametro di controllo della griglia di controllo (g1).
            Kg2 (float): Parametro di controllo della griglia schermo (g2).
            KvB (float): Parametro per il breakdown (solitamente per correnti di griglia).
        """
        if mu <= 0 or Kp <= 0 or X <= 0 or Kg1 <= 0 or Kg2 <= 0:
            raise ValueError("I parametri mu, Kp, X, Kg1, Kg2 devono essere positivi.")

        self.mu = float(mu)
        self.Kp = float(Kp)
        self.X = float(X)
        self.Kg1 = float(Kg1)
        self.Kg2 = float(Kg2) # Parametro specifico per la griglia schermo
        self.KvB = float(KvB)

        print(f"Pentodo creato con mu={self.mu}, Kp={self.Kp}, X={self.X}, Kg1={self.Kg1}, Kg2={self.Kg2}.")

    def calculate_anode_current(self, Vg1, Vg2, Va):
        """
        Calcola la corrente di Anodo (Ia) del Pentodo/Tetrodo.
        Estensione del modello di Koren che considera l'effetto di Vg2.
        Un approccio comune è usare una "tensione equivalente" che combina Vg1, Vg2 e Va.
        Questo modello è una semplificazione, i pentodi reali hanno un comportamento più complesso
        e richiederebbero modelli più robusti (es. di tipo "SPICE-like" o tabelle).
        """
        Ia = np.zeros_like(Vg1, dtype=float) if isinstance(Vg1, np.ndarray) else 0.0

        # Tensione equivalente per il controllo del flusso di elettroni (V_effective)
        # Questo è il cuore del modello del pentodo per come Vg1, Vg2 e Va influenzano Ia.
        # Formula di Koren per i pentodi usa spesso un V_equivalent combinato:
        # V_eq = (Va / mu_eff) + (Vg1 / Kg1) + (Vg2 / Kg2_eff)
        # dove mu_eff e Kg2_eff possono dipendere leggermente dalle tensioni.
        # Per un modello di base, possiamo semplificare:
        Ve = (Vg1 / self.Kg1) + (Vg2 / self.Kg2) # Griglie di controllo del flusso

        # Cutoff condition (quando Ve è sotto la soglia di conduzione)
        # La corrente è zero se la tensione effettiva è insufficiente.
        active_condition = Ve > 1e-6 # Assicuriamo che Ve sia positivo e non zero

        if np.any(active_condition):
            Ve_active = Ve[active_condition]
            # La corrente di anodo è calcolata come nel triodo, ma con Ve dipendente da Vg1 e Vg2
            Ia[active_condition] = self.Kp * (Ve_active ** self.X)

            # Aggiusta per saturazione del pentodo (Vds) o altre non linearità
            # I pentodi hanno una saturazione abbastanza piatta in regione attiva.
            # Questo modello semplice non cattura la "ginocchiatura" in regione di triodo.
            # Per una simulazione più fedele, ci vorrebbe un modello a più segmenti o un lookup table.
            # Qui si assume una saturazione ideale con Vg1 e Vg2 che dominano il controllo.

        # Assicurati che la corrente non sia mai negativa
        return np.maximum(0, Ia)

    def __str__(self):
        return (f"Pentode(mu={self.mu}, Kp={self.Kp}, X={self.X}, "
                f"Kg1={self.Kg1}, Kg2={self.Kg2})")

# Esempio di utilizzo e plot delle curve I-V
if __name__ == "__main__":
    # Parametri tipici per un EL34 o 6L6 (valvole di potenza)
    pentode1 = Pentode(mu=10.0, Kp=300.0, X=1.5, Kg1=5.0, Kg2=10.0)
    print(pentode1)

    import matplotlib.pyplot as plt

    # Plot delle curve Ia vs Va per diversi Vg1 (con Vg2 fisso)
    va_values = np.linspace(0, 400, 100) # Va da 0 a 400V
    vg1_steps = np.array([-5, -10, -15, -20]) # Diversi valori di Vg1 (griglia di controllo)
    vg2_fixed = 250.0 # Vg2 fissa (tensione della griglia schermo)

    plt.figure(figsize=(10, 6))
    for vg1 in vg1_steps:
        ia_values = pentode1.calculate_anode_current(vg1, vg2_fixed, va_values)
        plt.plot(va_values, ia_values * 1e3, label=f'Vg1={vg1:.1f}V') # Corrente in mA

    plt.title(f'Caratteristiche di Uscita del Pentodo (Ia vs Va) con Vg2={vg2_fixed}V')
    plt.xlabel('Va (V)')
    plt.ylabel('Ia (mA)')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot delle curve di trasferimento (Ia vs Vg1 per Va e Vg2 fissi)
    vg1_values_transfer = np.linspace(-30, 0, 100)
    va_fixed = 300.0 # Va fissa
    ia_transfer = pentode1.calculate_anode_current(vg1_values_transfer, vg2_fixed, va_fixed)

    plt.figure(figsize=(8, 5))
    plt.plot(vg1_values_transfer, ia_transfer * 1e3)
    plt.title(f'Caratteristiche di Trasferimento del Pentodo (Ia vs Vg1) con Va={va_fixed}V, Vg2={vg2_fixed}V')
    plt.xlabel('Vg1 (V)')
    plt.ylabel('Ia (mA)')
    plt.grid(True)
    plt.show()
