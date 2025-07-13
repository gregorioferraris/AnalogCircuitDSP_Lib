# components/capacitor.py

class Capacitor:
    def __init__(self, capacitance, sample_rate=48000):
        if capacitance <= 0:
            raise ValueError("La capacità deve essere un valore positivo.")
        if sample_rate <= 0:
            raise ValueError("La frequenza di campionamento deve essere un valore positivo.")

        self.C = float(capacitance)
        self.sample_rate = float(sample_rate)
        self.Ts = 1.0 / self.sample_rate # Periodo di campionamento

        # Stato interno per l'integrazione numerica (es. Trapezoidale)
        self.prev_voltage = 0.0
        self.prev_current = 0.0 # Potrebbe non essere usata direttamente, dipende dall'integrazione

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente del condensatore basandosi sulla variazione di tensione.
        Per un modello più accurato in un solutore di circuito, la corrente è
        implicita nella risoluzione del sistema di equazioni.
        Questo metodo è più per un calcolo diretto dV/dt in un contesto semplice.
        """
        # In un solutore MNA, la corrente è definita implicitamente dall'equazione
        # di integrazione trapezoidale che coinvolge lo stato precedente.
        # Questa implementazione è per un calcolo diretto della corrente data una tensione,
        # che potrebbe non essere direttamente usata nell'MNA per la corrente del condensatore.
        # Spesso, in MNA, si usa la forma: I_C(t) = C/Ts * (V_C(t) - V_C(t-Ts))
        # dove V_C(t-Ts) è prev_voltage.
        # Quindi, se usata, la corrente dipende anche dallo stato precedente.
        # Per semplicità, qui è la corrente in regime stazionario data una V, non una variazione.
        # L'MNA gestisce la dinamica.
        # Per ora restituisce 0, la logica dinamica è nel solutore MNA
        return 0.0 # La corrente dinamica viene gestita nel solutore tramite integrazione.


    def update_state(self, current_voltage, current_current=None):
        """
        Aggiorna lo stato interno del condensatore per il prossimo passo di simulazione.
        """
        self.prev_voltage = current_voltage
        if current_current is not None:
            self.prev_current = current_current

    def get_previous_voltage(self):
        """Restituisce la tensione ai capi del condensatore nel passo di tempo precedente."""
        return self.prev_voltage

    def __str__(self):
        return f"Capacitor({self.C*1e6:.2f} uF, Fs={self.sample_rate:.0f} Hz)"
