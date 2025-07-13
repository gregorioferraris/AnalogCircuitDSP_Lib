import numpy as np
from utils.constants import Q, K, T_ROOM, VT_ROOM # Importa costanti fisiche

class BJT:
    """
    Modello fisico semplificato per un BJT NPN (modello Ebers-Moll semplificato).
    Questo modello si concentra sulla relazione Ic-Vbe.
    """
    def __init__(self, Is=1e-14, Bf=100.0, Br=0.1, Vaf=np.inf, Var=np.inf, temperature_k=T_ROOM):
        """
        Inizializza il BJT.

        Args:
            Is (float): Corrente di saturazione del diodo di giunzione (A). Default 1e-14A.
            Bf (float): Guadagno di corrente in avanti (beta_forward). Default 100.
            Br (float): Guadagno di corrente inverso (beta_reverse). Default 0.1.
            Vaf (float): Tensione di Early in avanti (V). Default np.inf (nessun effetto Early).
            Var (float): Tensione di Early inversa (V). Default np.inf.
            temperature_k (float): Temperatura in Kelvin. Default temperatura ambiente.
        """
        if Is <= 0 or Bf <= 0 or Br <= 0 or temperature_k <= 0:
            raise ValueError("I parametri Is, Bf, Br, Temperature devono essere positivi.")

        self.Is = float(Is)
        self.Bf = float(Bf) # Beta Forward (hFE)
        self.Br = float(Br) # Beta Reverse
        self.Vaf = float(Vaf) # Forward Early Voltage
        self.Var = float(Var) # Reverse Early Voltage
        self.Vt = (K * temperature_k) / Q # Tensione termica

        # Correnti di saturazione inversa di diodi equivalenti per Ebers-Moll
        self.Ise = self.Is # Emitter-base saturation current
        self.Isc = self.Is # Collector-base saturation current

        print(f"BJT (NPN) creato con Is={self.Is}A, Bf={self.Bf}, Vt={self.Vt:.4f}V.")

    def calculate_collector_current(self, Vbe, Vce):
        """
        Calcola la corrente di Collettore (Ic) per un BJT NPN.
        (Modello semplificato, ignora la corrente di base e si concentra su Ic-Vbe)

        Args:
            Vbe (float or np.ndarray): Tensione Base-Emettitore (V).
            Vce (float or np.ndarray): Tensione Collettore-Emettitore (V).

        Returns:
            float or np.ndarray: Corrente di Collettore (Ic) (A).
        """
        # Equazione di Shockley per la giunzione BE
        # Assume che il BJT sia in regione attiva o saturazione.
        # La corrente di collettore è dominata dalla diffusione in avanti.
        if isinstance(Vbe, np.ndarray):
            Vbe_eff = np.clip(Vbe, None, 700 * self.Vt) # Protezione overflow exp
        else:
            Vbe_eff = min(Vbe, 700 * self.Vt)

        # Corrente di diffusione in avanti (relativa alla giunzione B-E)
        If = self.Is * (np.exp(Vbe_eff / self.Vt) - 1)

        # Effetto Early (modulazione della larghezza della base)
        # Aumenta la corrente di collettore all'aumentare di Vce
        early_factor = 1.0 + (Vce / self.Vaf) if self.Vaf != np.inf else 1.0

        # Corrente di collettore
        Ic = If * early_factor

        # Gestione delle regioni di Cutoff e Saturazione
        # In cutoff (Vbe basso), Ic è quasi 0. If già lo gestisce.
        # In saturazione, il modello semplificato non è preciso senza la corrente di base e di collettore-base inversa.
        # Per una simulazione più robusta, si userebbe un solutore di circuito completo
        # e un modello Ebers-Moll più completo con giunzioni CB e BE.
        # Per ora, limitiamo la corrente al minimo per evitare negative assurde.
        return np.maximum(0, Ic)

    def calculate_base_current(self, Vbe, Vce):
        """
        Calcola la corrente di Base (Ib) per un BJT NPN.
        (Semplificato: Ib = Ic / Bf)
        """
        Ic = self.calculate_collector_current(Vbe, Vce)
        return Ic / self.Bf

    def __str__(self):
        return f"BJT(Is={self.Is:.2e}A, Bf={self.Bf}, Vaf={self.Vaf}V)"

# Esempio di utilizzo e plot delle curve I-V
if __name__ == "__main__":
    bjt1 = BJT(Is=1e-14, Bf=200, Vaf=80.0) # BJT NPN tipico (es. 2N3904)
    print(bjt1)

    import matplotlib.pyplot as plt

    # Plot delle curve Ic vs Vce per diversi Vbe (o Ib)
    vce_values = np.linspace(0, 10, 100) # Vce da 0 a 10V
    vbe_steps = np.linspace(0.6, 0.75, 5) # Vbe da 0.6V a 0.75V (regione attiva)

    plt.figure(figsize=(10, 6))
    for vbe in vbe_steps:
        ic_values = bjt1.calculate_collector_current(vbe, vce_values)
        plt.plot(vce_values, ic_values * 1e3, label=f'Vbe={vbe:.2f}V') # Corrente in mA

    plt.title('Caratteristiche di Uscita del BJT (Ic vs Vce)')
    plt.xlabel('Vce (V)')
    plt.ylabel('Ic (mA)')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot delle curve di trasferimento (Ic vs Vbe per Vce fissa)
    vbe_values_transfer = np.linspace(0.5, 0.8, 100)
    vce_fixed = 5.0 # Vce fissa
    ic_transfer = bjt1.calculate_collector_current(vbe_values_transfer, vce_fixed)

    plt.figure(figsize=(8, 5))
    plt.plot(vbe_values_transfer, ic_transfer * 1e3)
    plt.title('Caratteristiche di Trasferimento del BJT (Ic vs Vbe)')
    plt.xlabel('Vbe (V)')
    plt.ylabel('Ic (mA)')
    plt.grid(True)
    plt.show()
