import numpy as np

class Triode:
    """
    Modello fisico per un Triodo (valvola termoionica a tre elettrodi)
    basato sul modello di Koren.
    """
    def __init__(self, mu=100.0, Kp=600.0, X=1.5, Kg1=2.5, KvB=300.0):
        """
        Inizializza il Triodo con i parametri del modello Koren.

        Args:
            mu (float): Fattore di amplificazione (mu).
            Kp (float): Parametro di transconduttanza (mA/V^X).
            X (float): Esponente non lineare.
            Kg1 (float): Parametro di controllo della griglia (influenza la conduttanza in avanti).
            KvB (float): Parametro per il breakdown in inversa (grid current, cathode emission).
                          Un valore alto (>100) significa che ignora l'effetto.
                          Per una simulazione pura dell'emissione termica, è più complesso.
        """
        if mu <= 0 or Kp <= 0 or X <= 0 or Kg1 <= 0:
            raise ValueError("I parametri mu, Kp, X, Kg1 devono essere positivi.")

        self.mu = float(mu)
        self.Kp = float(Kp)
        self.X = float(X)
        self.Kg1 = float(Kg1)
        self.KvB = float(KvB) # Per modellare la corrente di griglia (non primaria)

        print(f"Triodo creato con mu={self.mu}, Kp={self.Kp}, X={self.X}, Kg1={self.Kg1}.")

    def calculate_anode_current(self, Vg, Va):
        """
        Calcola la corrente di Anodo (Ia) del Triodo usando il modello di Koren.

        Args:
            Vg (float or np.ndarray): Tensione di Griglia rispetto al Catodo (Vgk) (V).
            Va (float or np.ndarray): Tensione di Anodo rispetto al Catodo (Vak) (V).

        Returns:
            float or np.ndarray: Corrente di Anodo (Ia) (A).
        """
        # Calcola la tensione equivalente (Ve) che controlla la corrente
        # Questa è la parte centrale del modello di Koren
        Ve = (Va / self.mu) + (Vg / self.Kg1)

        # Regola Ve per evitare valori negativi o troppo piccoli che non hanno senso
        # nel logaritmo o nella potenza frazionaria.
        # Una piccola offset può aiutare con la stabilità numerica.
        # La condizione Ve <= 0 è la regione di cutoff (o vicino)
        Ia = np.zeros_like(Vg, dtype=float) if isinstance(Vg, np.ndarray) else 0.0

        active_condition = Ve > 1e-6 # Assicuriamo che Ve sia positivo per evitare errori numerici

        if np.any(active_condition):
            Ve_active = Ve[active_condition]
            # La corrente di anodo è calcolata solo se Ve è maggiore di zero
            Ia[active_condition] = self.Kp * (Ve_active ** self.X)

        # La corrente di griglia non è inclusa esplicitamente qui,
        # ma in modelli più completi, la griglia può condurre (Grid Current)
        # se Vg diventa positiva rispetto al catodo (Vgk > 0.0).
        # Per ora, Ia sarà zero per Vg << 0 o Va << 0.

        # Assicurati che la corrente non sia mai negativa (non può fluire al contrario nell'anodo)
        return np.maximum(0, Ia)

    def __str__(self):
        return f"Triode(mu={self.mu}, Kp={self.Kp}, X={self.X}, Kg1={self.Kg1})"

# Esempio di utilizzo e plot delle curve I-V
if __name__ == "__main__":
    # Parametri tipici per un 12AX7 (un triodo molto comune)
    triode1 = Triode(mu=100.0, Kp=600.0, X=1.5, Kg1=2.5) # Corrisponde circa a una 12AX7
    print(triode1)

    import matplotlib.pyplot as plt

    # Plot delle curve Ia vs Va per diversi Vg
    va_values = np.linspace(0, 400, 100) # Va da 0 a 400V
    vg_steps = np.array([0, -0.5, -1, -2, -4, -8]) # Diversi valori di Vg (Vg=0 significa massima conduzione)

    plt.figure(figsize=(10, 6))
    for vg in vg_steps:
        ia_values = triode1.calculate_anode_current(vg, va_values)
        plt.plot(va_values, ia_values * 1e3, label=f'Vg={vg:.1f}V') # Corrente in mA

    plt.title('Caratteristiche di Uscita del Triodo (Ia vs Va)')
    plt.xlabel('Va (V)')
    plt.ylabel('Ia (mA)')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot delle curve di trasferimento (Ia vs Vg per Va fissa)
    vg_values_transfer = np.linspace(-10, 0, 100)
    va_fixed = 250.0 # Va fissa
    ia_transfer = triode1.calculate_anode_current(vg_values_transfer, va_fixed)

    plt.figure(figsize=(8, 5))
    plt.plot(vg_values_transfer, ia_transfer * 1e3)
    plt.title('Caratteristiche di Trasferimento del Triodo (Ia vs Vg)')
    plt.xlabel('Vg (V)')
    plt.ylabel('Ia (mA)')
    plt.grid(True)
    plt.show()
