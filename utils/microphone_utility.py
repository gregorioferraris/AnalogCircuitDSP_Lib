# utils/microphone_utility.py
import numpy as np
from scipy.signal import firwin, freqz, lfilter
import matplotlib.pyplot as plt

# Importa le classi dei microfoni
from components.condenser_microphone import CondenserMicrophone
from components.dynamic_microphone import DynamicMicrophone
from components.ribbon_microphone import RibbonMicrophone
from components.component import Component # Per il type hinting

class MicrophoneUtility:
    """
    Utility per calcolare la risposta in frequenza di un microfono
    e per esportarla come filtro FIR per l'applicazione in tempo reale.
    """
    def __init__(self, microphone: Component): # Usiamo Component per generalizzare
        """
        Inizializza l'utility con un'istanza di microfono.
        Args:
            microphone (Component): Un'istanza di CondenserMicrophone, DynamicMicrophone o RibbonMicrophone.
        """
        if not isinstance(microphone, (CondenserMicrophone, DynamicMicrophone, RibbonMicrophone)):
            raise TypeError("Il microfono deve essere un'istanza di CondenserMicrophone, DynamicMicrophone o RibbonMicrophone.")
        self.microphone = microphone
        self.fir_coeffs = None # Coefficienti FIR del filtro
        self.sample_rate = 44100 # Frequenza di campionamento di default, da impostare o passare

    def set_sample_rate(self, sample_rate: int):
        """Imposta la frequenza di campionamento per il calcolo del filtro."""
        self.sample_rate = sample_rate

    def calculate_frequency_response(self, frequencies: np.ndarray) -> np.ndarray:
        """
        Calcola la risposta in frequenza del microfono.
        Args:
            frequencies (np.ndarray): Array di frequenze in Hz.
        Returns:
            np.ndarray: Array di numeri complessi che rappresentano la risposta in frequenza (V/Pa).
        """
        # Ogni classe di microfono ha il suo metodo get_frequency_response
        return self.microphone.get_frequency_response(frequencies)

    def export_as_fir(self, frequency_response: np.ndarray, frequencies: np.ndarray, num_taps: int) -> np.ndarray:
        """
        Converte la risposta in frequenza in un filtro FIR (Finite Impulse Response).
        Args:
            frequency_response (np.ndarray): Risposta in frequenza complessa.
            frequencies (np.ndarray): Frequenze corrispondenti alla risposta (da 0 a Nyquist).
            num_taps (int): Numero di tap (lunghezza) del filtro FIR.
        Returns:
            np.ndarray: Coefficienti del filtro FIR.
        """
        if self.sample_rate is None or self.sample_rate <= 0:
            raise ValueError("Sample rate non impostato o non valido. Chiamare set_sample_rate prima.")

        nyquist_freq = self.sample_rate / 2.0
        
        # Interpolazione della risposta in frequenza per coprire l'intera banda di Nyquist
        n_fft = num_taps * 2 # Un numero di punti per l'IFFT, almeno 2*num_taps
        
        freq_interp = np.linspace(0, nyquist_freq, n_fft // 2 + 1)
        
        # Interpolazione della magnitudine (dB)
        magnitude_db = 20 * np.log10(np.abs(frequency_response + 1e-12)) # Aggiungi epsilon per log(0)
        interpolated_magnitude_db = np.interp(freq_interp, frequencies, magnitude_db)
        
        # Interpolazione della fase (radianti, unwrapped)
        phase_rad = np.unwrap(np.angle(frequency_response))
        interpolated_phase_rad = np.interp(freq_interp, frequencies, phase_rad)
        
        # Ricostruisci la risposta complessa interpolata
        interpolated_response_complex = 10**(interpolated_magnitude_db / 20) * np.exp(1j * interpolated_phase_rad)
        
        # Crea la risposta in frequenza completa per IFFT (simmetrica)
        full_spectrum = np.concatenate((interpolated_response_complex, np.conj(interpolated_response_complex[-2:0:-1])))
        
        # Calcola l'impulso di risposta (FIR coefficients)
        impulse_response = np.fft.ifft(full_spectrum)
        
        # Prendi la parte reale e tronca al numero di tap desiderato
        fir_coeffs = np.real(impulse_response[:num_taps])
        
        # Applica una finestra per ridurre gli artefatti (es. finestra di Hamming)
        window = np.hamming(num_taps)
        fir_coeffs = fir_coeffs * window
        
        # Normalizza i coefficienti per mantenere il guadagno (opzionale, dipende dall'uso)
        # fir_coeffs /= np.sum(fir_coeffs) # Solo se il guadagno DC deve essere 1

        self.fir_coeffs = fir_coeffs
        return fir_coeffs

    def process_audio_signal(self, input_acoustic_signal: np.ndarray) -> np.ndarray:
        """
        Applica il filtro FIR pre-calcolato a un segnale acustico in ingresso.
        Args:
            input_acoustic_signal (np.ndarray): Segnale di pressione acustica (es. da un file WAV).
        Returns:
            np.ndarray: Segnale elettrico in uscita dal microfono (Volt).
        """
        if self.fir_coeffs is None:
            raise ValueError("Il filtro FIR non è stato ancora calcolato. Chiamare export_as_fir prima.")
        
        # Applica il filtro FIR al segnale di input
        # lfilter è efficiente per l'applicazione di filtri FIR
        output_electrical_signal = lfilter(self.fir_coeffs, 1.0, input_acoustic_signal)
        return output_electrical_signal

    def plot_frequency_response(self, frequencies: np.ndarray, response: np.ndarray, title: str = "Risposta in Frequenza Microfono"):
        """
        Traccia la magnitudine e la fase della risposta in frequenza.
        Args:
            frequencies (np.ndarray): Frequenze in Hz.
            response (np.ndarray): Risposta in frequenza complessa.
            title (str): Titolo del grafico.
        """
        plt.figure(figsize=(12, 8))

        # Magnitudine
        plt.subplot(2, 1, 1)
        plt.semilogx(frequencies, 20 * np.log10(np.abs(response + 1e-12))) # Aggiungi epsilon per log(0)
        plt.title(f"{title} - Magnitudine")
        plt.xlabel("Frequenza (Hz)")
        plt.ylabel("Guadagno (dBV/Pa)") # Tensione per Pascal
        plt.grid(True, which="both", ls="-")
        plt.xlim([20, 20000]) # Limita l'asse X alle frequenze udibili
        plt.ylim([-80, 20]) # Range tipico per la sensibilità dei microfoni

        # Fase
        plt.subplot(2, 1, 2)
        plt.semilogx(frequencies, np.degrees(np.unwrap(np.angle(response))))
        plt.title(f"{title} - Fase")
        plt.xlabel("Frequenza (Hz)")
        plt.ylabel("Fase (gradi)")
        plt.grid(True, which="both", ls="-")
        plt.xlim([20, 20000])

        plt.tight_layout()
        plt.show()

# --- Esempio di utilizzo ---
if __name__ == "__main__":
    # Per questo esempio, assicurati che le classi dei microfoni siano importabili
    # o usa dei mock temporanei se non le hai ancora tutte.
    # Esempio di classi microfono minimali per il test
    class MockCondenserMicrophone(Component):
        def __init__(self, name, out_node, gnd_node):
            super().__init__(name, out_node, gnd_node)
        def get_frequency_response(self, frequencies):
            # Risposta piatta con un roll-off in bassa frequenza
            omega = 2 * np.pi * frequencies
            Rc = 1e9 # Impedenza di ingresso preamp
            Cc = 50e-12 # Capacità capsula
            sensitivity = 0.01 # 10 mV/Pa
            return sensitivity * (1j * omega * Rc * Cc) / (1 + 1j * omega * Rc * Cc)
    class MockDynamicMicrophone(Component):
        def __init__(self, name, out_node, gnd_node):
            super().__init__(name, out_node, gnd_node)
        def get_frequency_response(self, frequencies):
            # Risposta con risonanza in media frequenza e roll-off in HF/LF
            omega = 2 * np.pi * frequencies
            f_res = 1000 # Hz
            Q = 1.5
            sensitivity = 0.002 # 2 mV/Pa
            # Modello di risonatore di secondo ordine
            response = sensitivity / (1 + 1j * Q * (frequencies/f_res - f_res/frequencies))
            return response
    class MockRibbonMicrophone(Component):
        def __init__(self, name, out_node, gnd_node):
            super().__init__(name, out_node, gnd_node)
        def get_frequency_response(self, frequencies):
            # Risposta con roll-off in alta frequenza (trasformatore) e risonanza bassa
            omega = 2 * np.pi * frequencies
            f_res = 50 # Hz
            Q = 0.8
            sensitivity = 0.0015 # 1.5 mV/Pa
            # Risonanza meccanica
            mech_response = 1.0 / (1 + 1j * Q * (frequencies/f_res - f_res/frequencies))
            # Roll-off in alta frequenza (es. filtro passa-basso di primo ordine)
            f_lp = 10000 # Hz
            lp_response = 1.0 / (1 + 1j * (frequencies / f_lp))
            return sensitivity * mech_response * lp_response

    # Sostituisci gli import reali con i mock per questo test se non hai ancora tutti i file reali
    globals().update({
        "CondenserMicrophone": MockCondenserMicrophone,
        "DynamicMicrophone": MockDynamicMicrophone,
        "RibbonMicrophone": MockRibbonMicrophone
    })

    sample_rate = 44100 # Hz
    num_fir_taps = 512 # Lunghezza del filtro FIR

    # --- Esempio 1: Microfono a Condensatore ---
    condenser_mic = CondenserMicrophone("MyCondenser", "mic_out", "GND")
    mic_utility_condenser = MicrophoneUtility(condenser_mic)
    mic_utility_condenser.set_sample_rate(sample_rate)

    frequencies = np.logspace(np.log10(20), np.log10(20000), 500)
    response_condenser = mic_utility_condenser.calculate_frequency_response(frequencies)
    mic_utility_condenser.plot_frequency_response(frequencies, response_condenser, "Risposta Microfono a Condensatore")

    fir_coeffs_condenser = mic_utility_condenser.export_as_fir(response_condenser, frequencies, num_fir_taps)
    print(f"\nCondenser Mic FIR Coeffs ({len(fir_coeffs_condenser)} tap): {fir_coeffs_condenser[:5]}...")

    # --- Esempio 2: Microfono Dinamico ---
    dynamic_mic = DynamicMicrophone("MyDynamic", "mic_out", "GND")
    mic_utility_dynamic = MicrophoneUtility(dynamic_mic)
    mic_utility_dynamic.set_sample_rate(sample_rate)

    response_dynamic = mic_utility_dynamic.calculate_frequency_response(frequencies)
    mic_utility_dynamic.plot_frequency_response(frequencies, response_dynamic, "Risposta Microfono Dinamico")

    fir_coeffs_dynamic = mic_utility_dynamic.export_as_fir(response_dynamic, frequencies, num_fir_taps)
    print(f"\nDynamic Mic FIR Coeffs ({len(fir_coeffs_dynamic)} tap): {fir_coeffs_dynamic[:5]}...")

    # --- Esempio 3: Microfono a Nastro ---
    ribbon_mic = RibbonMicrophone("MyRibbon", "mic_out", "GND")
    mic_utility_ribbon = MicrophoneUtility(ribbon_mic)
    mic_utility_ribbon.set_sample_rate(sample_rate)

    response_ribbon = mic_utility_ribbon.calculate_frequency_response(frequencies)
    mic_utility_ribbon.plot_frequency_response(frequencies, response_ribbon, "Risposta Microfono a Nastro")

    fir_coeffs_ribbon = mic_utility_ribbon.export_as_fir(response_ribbon, frequencies, num_fir_taps)
    print(f"\nRibbon Mic FIR Coeffs ({len(fir_coeffs_ribbon)} tap): {fir_coeffs_ribbon[:5]}...")

    # --- Come usare in un loop di simulazione (concettuale) ---
    print("\n--- Esempio di integrazione in un loop di simulazione ---")
    # Immagina di avere un segnale di pressione acustica (es. da un'onda sonora simulata)
    acoustic_pressure_signal = np.random.randn(sample_rate) * 0.5 # 1 secondo di "pressione acustica"

    # Processa il segnale acustico con l'utility del microfono
    electrical_output_signal = mic_utility_condenser.process_audio_signal(acoustic_pressure_signal)

    print(f"Segnale elettrico di output del microfono (primi 5 campioni): {electrical_output_signal[:5]}...")

    # Questo 'electrical_output_signal' può essere usato per aggiornare una VoltageSource
    # nel tuo circuito MNA ad ogni passo temporale.
    # Esempio concettuale nel tuo MnaSolver.solve_transient loop:
    # for i, t in enumerate(time_points):
    #     current_acoustic_pressure = acoustic_pressure_signal[i] # Preleva il campione di pressione
    #     mic_electrical_output = mic_utility_condenser.process_audio_signal(np.array([current_acoustic_pressure]))[0]
    #     # Trova la VoltageSource nel circuito che rappresenta l'uscita del microfono
    #     # mic_output_vs = next((c for c in self.circuit.get_voltage_sources() if c.name == "MicOutputVS"), None)
    #     # if mic_output_vs:
    #     #     mic_output_vs.set_voltage(mic_electrical_output)
    #     # ... poi continua con il normale processo di risoluzione MNA ...

