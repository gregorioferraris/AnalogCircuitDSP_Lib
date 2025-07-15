# components/tape.py
import numpy as np
from components.component import Component

class Tape(Component):
    def __init__(self, name: str, input_node: str, output_node: str, ground_node: str,
                 saturation_level: float = 1.0, # Livello di saturazione (es. in Weber/metro)
                 gain: float = 1.0,             # Guadagno lineare complessivo del nastro
                 bias_frequency: float = 150e3, # Frequenza di bias (per linearizzazione e riduzione rumore)
                 tape_speed: float = 0.38,      # Velocità del nastro (m/s) (es. 0.38 m/s per 15 IPS)
                 noise_level: float = 0.0,      # Livello di rumore additivo (RMS)
                 hysteresis_factor: float = 0.0 # Fattore per modellare l'isteresi (semplificato)
                 ):
        """
        Inizializza un modello semplificato di nastro magnetico.
        Questo modello mira a catturare effetti come la saturazione.
        Args:
            name (str): Nome univoco dell'istanza (es. "TAPE1").
            input_node (str): Nome del nodo di ingresso del segnale.
            output_node (str): Nome del nodo di uscita del segnale.
            ground_node (str): Nome del nodo di riferimento (ground).
            saturation_level (float): Livello a cui il nastro inizia a saturare (es. in unità arbitrarie di flusso magnetico).
            gain (float): Guadagno lineare applicato prima della non linearità.
            bias_frequency (float): Frequenza del segnale di bias AC, usato per linearizzare il nastro.
            tape_speed (float): Velocità lineare del nastro in metri al secondo (m/s).
            noise_level (float): Livello RMS di rumore bianco additivo (da aggiungere all'uscita).
            hysteresis_factor (float): Un fattore per modellare l'isteresi (0 per nessuno, >0 per isteresi).
        """
        super().__init__(name, input_node, output_node, ground_node)
        self.pin_names = ('input', 'output', 'ground')

        self.saturation_level = saturation_level
        self.gain = gain
        self.bias_frequency = bias_frequency
        self.tape_speed = tape_speed
        self.noise_level = noise_level
        self.hysteresis_factor = hysteresis_factor

        # Stato interno per l'isteresi o per effetti dinamici
        self._prev_input_signal = 0.0
        self._prev_output_signal = 0.0

    def calculate_output(self, input_signal: float) -> float:
        """
        Calcola il segnale di uscita del nastro, includendo saturazione e una forma semplificata di isteresi.
        Questo modello è comportamentale e non basato su nodi MNA tradizionali.
        Args:
            input_signal (float): Segnale di tensione in ingresso al nastro.
        Returns:
            float: Segnale di tensione in uscita dal nastro.
        """
        # Applica guadagno lineare
        linear_signal = input_signal * self.gain

        # Modello di saturazione (funzione tangente iperbolica per una curva smooth)
        # La tanh è una buona approssimazione della curva di magnetizzazione
        saturated_signal = self.saturation_level * np.tanh(linear_signal / self.saturation_level)

        # Semplice modello di isteresi (dipende dal segnale precedente e dalla direzione del cambiamento)
        # Questo è un modello euristico, non una modellazione fisica completa dell'isteresi magnetica.
        if self.hysteresis_factor > 0:
            delta_input = input_signal - self._prev_input_signal
            if delta_input > 0: # Segnale in aumento
                saturated_signal = saturated_signal + self.hysteresis_factor * (saturated_signal - self._prev_output_signal)
            elif delta_input < 0: # Segnale in diminuzione
                saturated_signal = saturated_signal - self.hysteresis_factor * (self._prev_output_signal - saturated_signal)
            # Clip per mantenere entro i limiti di saturazione
            saturated_signal = np.clip(saturated_signal, -self.saturation_level, self.saturation_level)

        # Aggiunta di rumore (semplice rumore bianco)
        noise = np.random.normal(0, self.noise_level)
        output_signal = saturated_signal + noise

        # Aggiorna lo stato per il prossimo passo
        self._prev_input_signal = input_signal
        self._prev_output_signal = output_signal
        
        return output_signal

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Il nastro magnetico è un blocco funzionale. La sua uscita è calcolata comportamentalmente
        e dovrà essere applicata come sorgente di tensione controllata nel vettore RHS dal solutore.
        Per il metodo MNA, restituisce stamp_A e stamp_B vuoti, poiché il suo contributo diretto
        non è una semplice relazione lineare G*V = I.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

