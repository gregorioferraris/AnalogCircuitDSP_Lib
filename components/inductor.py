import numpy as np

class Inductor:
    """
    Modello fisico per un induttore usando il metodo trapezoidale (bilinear transform).
    """
    def __init__(self, inductance=10e-3, sample_rate=44100):
        """
        Inizializza l'induttore.

        Args:
            inductance (float): Il valore di induttanza in Henry. Default 10 mH.
            sample_rate (float): La frequenza di campionamento in Hz.
        """
        if inductance <= 0:
            raise ValueError("L'induttanza deve essere un valore positivo.")
        if sample_rate <= 0:
            raise ValueError("La frequenza di campionamento deve essere positiva.")

        self.L = float(inductance)
        self.Ts = 1.0 / sample_rate # Periodo di campionamento

        # Variabili di stato per il metodo trapezoidale
        # Voltage = (2L/Ts) * (I_now - I_prev) - V_prev
        # O, in termini di resistenza equivalente per la corrente: V_now = R_eq * I_now - V_prev_eff
        self.R_eq = (2.0 * self.L) / self.Ts # Resistenza equivalente
        self.current_prev = 0.0 # Corrente attraverso l'induttore al passo precedente
        self.voltage_prev = 0.0 # Tensione ai capi dell'induttore al passo precedente

        print(f"Induttore creato con L = {self.L}H, Fs = {sample_rate}Hz.")

    def calculate_current(self, voltage_now):
        """
        Calcola la corrente che scorre attraverso l'induttore al tempo attuale,
        data la tensione ai suoi capi.
        Questa funzione è per l'integrazione in un solutore di circuito nodale.
        La tensione_now è quella che stiamo cercando di risolvere.

        Args:
            voltage_now (float or np.ndarray): Tensione (V) ai capi dell'induttore al campione attuale.

        Returns:
            float or np.ndarray: Corrente (A) attraverso l'induttore.
        """
        # Equazione discreta dell'induttore (metodo trapezoidale, riarrangiata per corrente):
        # I_n = (Ts / (2L)) * (V_n + V_{n-1}) + I_{n-1}
        # In questo caso, per l'analisi nodale, è più utile esprimere la corrente
        # in funzione della tensione corrente e dello stato passato.
        # I_now = (1/R_eq) * (voltage_now - (-voltage_prev_eff))
        # dove -voltage_prev_eff = - (voltage_prev + self.R_eq * self.current_prev)
        current_now = (self.Ts / (2.0 * self.L)) * (voltage_now + self.voltage_prev) + self.current_prev
        return current_now


    def update_state(self, voltage_now, current_now):
        """
        Aggiorna le variabili di stato (tensione e corrente precedenti) per il prossimo passo.
        Deve essere chiamato dopo che il solutore ha determinato la tensione e la corrente del campione attuale.

        Args:
            voltage_now (float): Tensione (V) ai capi dell'induttore al campione attuale.
            current_now (float): Corrente (A) attraverso l'induttore al campione attuale.
        """
        self.voltage_prev = voltage_now
        self.current_prev = current_now

    def __str__(self):
        return f"Induttore(L={self.L}H, Fs={1.0/self.Ts}Hz)"

# Esempio di utilizzo (per un test isolato, non in un circuito completo)
if __name__ == "__main__":
    fs = 48000 # Frequenza di campionamento
    l1 = Inductor(10e-3, sample_rate=fs) # 10 mH

    # Simuliamo un segnale sinusoidale di tensione per vedere la corrente
    time_steps = np.arange(0, 0.01, l1.Ts) # 10ms di simulazione
    input_voltage = 1.0 * np.sin(2 * np.pi * 500 * time_steps) # 500Hz sinusoidale, 1Vpk

    output_current = []
    for i, v_now in enumerate(input_voltage):
        # In un vero solutore, qui risolveremmo il circuito per v_now e i_now
        # Per questo esempio, assumiamo v_now è nota e calcoliamo i_now direttamente
        i_now = l1.calculate_current(v_now)
        output_current.append(i_now)
        l1.update_state(v_now, i_now) # Aggiorniamo lo stato per il prossimo passo

    print(f"Prime 5 correnti calcolate: {output_current[:5]}")
