import numpy as np

class Capacitor:
    """
    Modello fisico per un condensatore usando il metodo trapezoidale (bilinear transform).
    """
    def __init__(self, capacitance=100e-9, sample_rate=44100):
        """
        Inizializza il condensatore.

        Args:
            capacitance (float): Il valore di capacità in Farad. Default 100 nF.
            sample_rate (float): La frequenza di campionamento in Hz.
        """
        if capacitance <= 0:
            raise ValueError("La capacità deve essere un valore positivo.")
        if sample_rate <= 0:
            raise ValueError("La frequenza di campionamento deve essere positiva.")

        self.C = float(capacitance)
        self.Ts = 1.0 / sample_rate  # Periodo di campionamento

        # Variabili di stato per il metodo trapezoidale
        # Current = (2C/Ts) * (V_now - V_prev) - I_prev
        # O, in termini di conduttanza equivalente: I_now = G_eq * V_now - I_prev_eff
        self.G_eq = (2.0 * self.C) / self.Ts # Conduttanza equivalente
        self.voltage_prev = 0.0 # Tensione ai capi del condensatore al passo precedente
        self.current_prev = 0.0 # Corrente attraverso il condensatore al passo precedente

        print(f"Condensatore creato con C = {self.C}F, Fs = {sample_rate}Hz.")

    def calculate_current(self, voltage_now):
        """
        Calcola la corrente che scorre attraverso il condensatore al tempo attuale,
        data la tensione ai suoi capi.
        Questa funzione è per l'integrazione in un solutore di circuito nodale.
        La tensione_now è quella che stiamo cercando di risolvere.

        Args:
            voltage_now (float or np.ndarray): Tensione (V) ai capi del condensatore al campione attuale.

        Returns:
            float or np.ndarray: Corrente (A) attraverso il condensatore.
        """
        # Equazione discreta del condensatore (metodo trapezoidale):
        # I_n = (2C/Ts) * (V_n - V_{n-1}) - I_{n-1}
        # Dove V_n è voltage_now, V_{n-1} è self.voltage_prev, I_{n-1} è self.current_prev
        current_now = self.G_eq * (voltage_now - self.voltage_prev) - self.current_prev
        return current_now

    def update_state(self, voltage_now, current_now):
        """
        Aggiorna le variabili di stato (tensione e corrente precedenti) per il prossimo passo.
        Deve essere chiamato dopo che il solutore ha determinato la tensione e la corrente del campione attuale.

        Args:
            voltage_now (float): Tensione (V) ai capi del condensatore al campione attuale.
            current_now (float): Corrente (A) attraverso il condensatore al campione attuale.
        """
        self.voltage_prev = voltage_now
        self.current_prev = current_now

    def __str__(self):
        return f"Condensatore(C={self.C}F, Fs={1.0/self.Ts}Hz)"

# Esempio di utilizzo (per un test isolato, non in un circuito completo)
if __name__ == "__main__":
    fs = 48000 # Frequenza di campionamento
    c1 = Capacitor(10e-9, sample_rate=fs) # 10 nF

    # Simuliamo un segnale sinusoidale di tensione per vedere la corrente
    # In un circuito reale, voltage_now verrebbe dal solutore.
    # Qui usiamo un esempio semplificato per mostrare l'aggiornamento dello stato.
    time_steps = np.arange(0, 0.01, c1.Ts) # 10ms di simulazione
    input_voltage = 1.0 * np.sin(2 * np.pi * 1000 * time_steps) # 1kHz sinusoidale, 1Vpk

    output_current = []
    for i, v_now in enumerate(input_voltage):
        # In un vero solutore, qui risolveremmo il circuito per v_now e i_now
        # Per questo esempio, assumiamo v_now è nota e calcoliamo i_now direttamente
        i_now = c1.calculate_current(v_now)
        output_current.append(i_now)
        c1.update_state(v_now, i_now) # Aggiorniamo lo stato per il prossimo passo

    print(f"Prime 5 correnti calcolate: {output_current[:5]}")
    # Nota: per vedere un comportamento significativo, avresti bisogno di un solutore di circuito completo.
