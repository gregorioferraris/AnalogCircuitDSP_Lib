# components/inductor.py

class Inductor:
    def __init__(self, inductance, sample_rate=48000):
        if inductance <= 0:
            raise ValueError("L'induttanza deve essere un valore positivo.")
        if sample_rate <= 0:
            raise ValueError("La frequenza di campionamento deve essere un valore positivo.")

        self.L = float(inductance)
        self.sample_rate = float(sample_rate)
        self.Ts = 1.0 / self.sample_rate # Periodo di campionamento

        # Stato interno per l'integrazione numerica (es. Trapezoidale)
        self.prev_current = 0.0
        self.prev_voltage = 0.0 # Potrebbe non essere usata direttamente, dipende dall'integrazione

    def calculate_current(self, voltage_difference):
        """
        Calcola la corrente che scorre attraverso l'induttore basandosi sulla tensione.
        Per un modello più accurato in un solutore di circuito, la corrente è
        implicita nella risoluzione del sistema di equazioni differenziali.
        """
        # In un solutore MNA, la corrente è una variabile ausiliaria e la sua dinamica
        # è gestita implicitamente tramite l'equazione dell'induttore e lo stato precedente.
        # Questa implementazione è per un calcolo diretto V = L di/dt, non direttamente usato nell'MNA.
        # Per ora restituisce 0, la logica dinamica è nel solutore MNA
        return 0.0 # La corrente dinamica viene gestita nel solutore tramite integrazione.

    def update_state(self, current_voltage=None, current_current=None):
        """
        Aggiorna lo stato interno dell'induttore per il prossimo passo di simulazione.
        """
        if current_current is not None:
            self.prev_current = current_current
        if current_voltage is not None:
            self.prev_voltage = current_voltage

    def get_previous_current(self):
        """Restituisce la corrente attraverso l'induttore nel passo di tempo precedente."""
        return self.prev_current

    def __str__(self):
        return f"Inductor({self.L*1e3:.2f} mH, Fs={self.sample_rate:.0f} Hz)"
