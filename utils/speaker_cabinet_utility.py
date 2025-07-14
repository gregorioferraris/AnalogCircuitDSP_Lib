# utils/speaker_cabinet_utility.py
import numpy as np
from scipy.signal import firwin, freqz, lfilter
import matplotlib.pyplot as plt

# Importa le classi dei componenti che useremo
from components.speaker_driver import SpeakerDriver
from components.closed_box_cabinet import ClosedBoxCabinet
from components.bass_reflex_cabinet import BassReflexCabinet

class SpeakerCabinetUtility:
    """
    Utility per calcolare la risposta in frequenza combinata di un driver e un cabinet,
    e per esportarla come filtro FIR per l'applicazione in tempo reale.
    """
    def __init__(self, driver: SpeakerDriver, cabinet: [ClosedBoxCabinet, BassReflexCabinet]):
        """
        Inizializza l'utility con un driver e un cabinet.
        Args:
            driver (SpeakerDriver): L'istanza del driver dell'altoparlante.
            cabinet (Union[ClosedBoxCabinet, BassReflexCabinet]): L'istanza del cabinet.
        """
        self.driver = driver
        self.cabinet = cabinet

        # Verifica la compatibilità dei nodi
        if self.driver.node_names_str[2] != self.cabinet.node_names_str[0]:
            print(f"Warning: Il nodo di velocità meccanica del driver ('{self.driver.node_names_str[2]}') "
                  f"non corrisponde al nodo di input del cabinet ('{self.cabinet.node_names_str[0]}'). "
                  f"Assicurati che siano collegati correttamente nel modello concettuale.")

    def calculate_combined_frequency_response(self, frequencies: np.ndarray, driver_input_voltage: float = 1.0) -> np.ndarray:
        """
        Calcola la risposta in frequenza combinata del sistema driver-cabinet.
        Assumiamo che l'output sia la pressione sonora generata dalla velocità del cono.
        
        La pressione sonora è proporzionale al flusso volumetrico (velocità del cono * area del cono).
        Quindi, calcoliamo la velocità del cono per una data tensione in ingresso.
        
        Args:
            frequencies (np.ndarray): Array di frequenze in Hz su cui calcolare la risposta.
            driver_input_voltage (float): Tensione RMS di ingresso al driver (per normalizzazione).
        Returns:
            np.ndarray: Array di numeri complessi che rappresentano la risposta in frequenza (guadagno e fase).
        """
        combined_response = np.zeros_like(frequencies, dtype=complex)
        
        for i, f in enumerate(frequencies):
            if f == 0: # Gestione DC
                # A DC, gli induttori sono cortocircuiti, i condensatori sono circuiti aperti.
                # La velocità del cono a DC è 0 (non c'è movimento continuo).
                combined_response[i] = 0.0
                continue

            # 1. Impedenza elettrica del driver (Re + jwLe)
            Z_elec = self.driver.get_electrical_impedance(f)

            # 2. Impedenza acustica del cabinet (analogia elettrica)
            Z_acoustic_cabinet = self.cabinet.get_acoustic_impedance_analogy(f)

            # 3. Impedenza meccanica del driver (in analogia elettrica)
            Z_mech_driver = self.driver.get_mechanical_impedance_analogy(f)

            # 4. Combinazione dell'impedenza meccanica del driver e del cabinet
            # Il driver "vede" l'impedenza del cabinet in parallelo con la sua propria impedenza di radiazione acustica
            # e l'impedenza del cabinet.
            # In questa analogia, il driver e il cabinet sono in parallelo dal punto di vista meccanico/acustico.
            # 1/Z_total_mech = 1/Z_mech_driver + 1/Z_acoustic_cabinet
            Y_total_mech = (1.0 / Z_mech_driver) + (1.0 / Z_acoustic_cabinet)
            Z_total_mech = 1.0 / Y_total_mech if Y_total_mech != 0 else np.inf

            # 5. Calcolo della funzione di trasferimento dalla tensione elettrica alla velocità del cono
            # V_cone_velocity / V_electrical_input = Bl / ( (Re + jwLe) * Z_total_mech + Bl^2 )
            # Questo è il guadagno dalla tensione elettrica alla velocità meccanica (tensione analogica)
            denominator = (Z_elec * Z_total_mech) + (self.driver.Bl**2)
            if denominator == 0:
                combined_response[i] = np.inf # Risonanza infinita
            else:
                velocity_gain = self.driver.Bl / denominator
                
                # La pressione sonora è proporzionale alla velocità del cono * area del cono
                # P = rho * c * V_cone_velocity * Sd (per un'onda piana semplice)
                # Per un altoparlante, è più complesso, ma la velocità del cono è il proxy.
                # Qui, normalizziamo la velocità per la tensione di ingresso e l'area del cono per ottenere un "guadagno acustico"
                acoustic_output_gain = velocity_gain * self.driver.Sd # Velocità del cono * Area -> Flusso volumetrico
                
                combined_response[i] = acoustic_output_gain
                
        return combined_response / driver_input_voltage # Normalizza per la tensione di ingresso

    def export_as_fir(self, frequency_response: np.ndarray, frequencies: np.ndarray, sample_rate: int, num_taps: int) -> np.ndarray:
        """
        Converte la risposta in frequenza in un filtro FIR (Finite Impulse Response).
        Args:
            frequency_response (np.ndarray): Risposta in frequenza complessa (dal metodo precedente).
            frequencies (np.ndarray): Frequenze corrispondenti alla risposta (da 0 a Nyquist).
            sample_rate (int): Frequenza di campionamento per il filtro FIR.
            num_taps (int): Numero di tap (lunghezza) del filtro FIR. Più tap = più fedeltà, più computazione.
        Returns:
            np.ndarray: Coefficienti del filtro FIR.
        """
        # Creiamo un array di frequenze normalizzate (da 0 a 1, dove 1 è la frequenza di Nyquist)
        nyquist_freq = sample_rate / 2.0
        normalized_frequencies = frequencies / nyquist_freq

        # Interpoliamo la risposta in frequenza per coprire l'intera banda di Nyquist
        # (se frequencies non copre l'intervallo 0-Nyquist in modo uniforme)
        # Per semplicità, assumiamo che frequencies sia già ben distribuito.
        
        # Usiamo il metodo "window" per la progettazione FIR (es. firwin2 o remez)
        # scipy.signal.firwin2 è più adatto per risposte arbitrarie.
        # Tuttavia, per una risposta in frequenza complessa, è più comune usare la trasformata inversa di Fourier.
        
        # Per la conversione da risposta in frequenza a FIR:
        # 1. Estendi la risposta in frequenza per essere simmetrica per una trasformata reale.
        # 2. Esegui una IFFT (Inverse Fast Fourier Transform).
        # 3. Applica un finestratura (windowing) e tronca per ottenere il numero desiderato di tap.

        # Assicurati che la risposta in frequenza sia definita fino alla frequenza di Nyquist.
        # Se frequencies non include 0 Hz e Nyquist, dovremmo aggiungerli.
        
        # Costruisci la risposta in frequenza completa per IFFT (simmetrica)
        # Frequenze da 0 a Fs/2 (Nyquist)
        # Numero di punti FFT = (num_taps - 1) * 2 (per simmetria) o una potenza di 2
        
        # Per semplicità, useremo un approccio che si basa su fft.ifft
        # Interpoliamo la risposta per avere un numero di punti adeguato per l'IFFT
        n_fft = num_taps * 2 # Un numero di punti per l'IFFT, almeno 2*num_taps
        
        # Crea un array di frequenze uniformemente spaziate per l'IFFT
        freq_interp = np.linspace(0, nyquist_freq, n_fft // 2 + 1)
        
        # Interpolazione della risposta in frequenza
        # np.interp funziona solo con valori reali. Dobbiamo interpolare magnitudine e fase separatamente
        # o usare un'interpolazione complessa.
        # Per semplicità, useremo np.interp per la magnitudine e np.unwrap per la fase.
        
        # Interpolazione della magnitudine (dB)
        magnitude_db = 20 * np.log10(np.abs(frequency_response))
        interpolated_magnitude_db = np.interp(freq_interp, frequencies, magnitude_db)
        
        # Interpolazione della fase (radianti, unwrapped)
        phase_rad = np.unwrap(np.angle(frequency_response))
        interpolated_phase_rad = np.interp(freq_interp, frequencies, phase_rad)
        
        # Ricostruisci la risposta complessa interpolata
        interpolated_response_complex = 10**(interpolated_magnitude_db / 20) * np.exp(1j * interpolated_phase_rad)
        
        # Crea la risposta in frequenza completa per IFFT (simmetrica)
        # H[k] = H*[N-k] per segnali reali
        full_spectrum = np.concatenate((interpolated_response_complex, np.conj(interpolated_response_complex[-2:0:-1])))
        
        # Calcola l'impulso di risposta (FIR coefficients)
        impulse_response = np.fft.ifft(full_spectrum)
        
        # Prendi la parte reale e tronca al numero di tap desiderato
        fir_coeffs = np.real(impulse_response[:num_taps])
        
        # Applica una finestra per ridurre gli artefatti (es. finestra di Hamming)
        window = np.hamming(num_taps)
        fir_coeffs = fir_coeffs * window
        
        # Normalizza i coefficienti per mantenere il guadagno (opzionale, dipende dall'uso)
        # fir_coeffs /= np.sum(fir_coeffs)

        return fir_coeffs

    def plot_frequency_response(self, frequencies: np.ndarray, response: np.ndarray, title: str = "Risposta in Frequenza"):
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
        plt.semilogx(frequencies, 20 * np.log10(np.abs(response)))
        plt.title(f"{title} - Magnitudine")
        plt.xlabel("Frequenza (Hz)")
        plt.ylabel("Guadagno (dB)")
        plt.grid(True, which="both", ls="-")
        plt.xlim([20, 20000]) # Limita l'asse X alle frequenze udibili

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
    # Parametri di esempio per un driver e un cabinet
    # Parametri Thiele-Small tipici per un woofer da 12 pollici
    # Re = 6.3 Ohm, Le = 0.5 mH, Bl = 12.0 Tm, Mms = 0.04 kg, Cms = 0.0003 m/N, Rms = 1.5 Ns/m, Sd = 0.05 m^2
    driver = SpeakerDriver(
        name="MyWoofer",
        elec_plus_node="elec_in_plus",
        elec_minus_node="elec_in_minus",
        mech_velocity_node="mech_vel",
        Re=6.3, Le=0.5e-3, Bl=12.0, Mms=0.04, Cms=0.0003, Rms=1.5, Sd=0.05
    )

    # Cabinet a cassa chiusa (Closed Box)
    # Volume tipico per un 12 pollici: 40 litri = 0.04 m^3
    closed_cabinet = ClosedBoxCabinet(
        name="MyClosedBox",
        mech_velocity_input_node="mech_vel",
        box_volume=0.04
    )

    # Cabinet Bass Reflex
    # Volume: 40 litri, Area port: 0.005 m^2 (es. 8cm diametro), Lunghezza port: 0.1 m
    bass_reflex_cabinet = BassReflexCabinet(
        name="MyBassReflex",
        mech_velocity_input_node="mech_vel",
        box_volume=0.04,
        port_area=0.005,
        port_length=0.1
    )

    # Inizializza l'utility con il driver e un cabinet (es. cassa chiusa)
    speaker_system = SpeakerCabinetUtility(driver, closed_cabinet)
    # speaker_system = SpeakerCabinetUtility(driver, bass_reflex_cabinet) # Per testare il bass reflex

    # Frequenze per il calcolo della risposta (da 20 Hz a 20 kHz, logaritmico)
    frequencies = np.logspace(np.log10(20), np.log10(20000), 500)

    # Calcola la risposta in frequenza combinata
    response = speaker_system.calculate_combined_frequency_response(frequencies)

    # Traccia la risposta
    speaker_system.plot_frequency_response(frequencies, response, "Risposta Driver + Closed Box")

    # Esempio di esportazione come filtro FIR
    sample_rate = 44100 # Hz
    num_fir_taps = 256 # Lunghezza del filtro FIR

    # Per l'esportazione FIR, è meglio avere frequenze fino a Nyquist
    fir_frequencies = np.linspace(0, sample_rate / 2, 500)
    fir_response = speaker_system.calculate_combined_frequency_response(fir_frequencies)
    fir_coeffs = speaker_system.export_as_fir(fir_response, fir_frequencies, sample_rate, num_fir_taps)

    print(f"\nCoefficienti FIR generati ({len(fir_coeffs)} tap):")
    print(fir_coeffs[:10], "...") # Stampa solo i primi 10 coefficienti

    # --- Esempio di applicazione del filtro FIR a un segnale audio (simulato) ---
    # Genera un segnale di rumore bianco di test
    test_signal_length = sample_rate # 1 secondo di rumore bianco
    test_signal = np.random.randn(test_signal_length) * 0.1

    # Applica il filtro FIR
    filtered_signal = lfilter(fir_coeffs, 1.0, test_signal)

    print(f"\nSegnale filtrato (primi 10 campioni): {filtered_signal[:10]}")

    # Plot del segnale filtrato (opzionale)
    # plt.figure(figsize=(10, 4))
    # plt.plot(np.arange(len(filtered_signal)) / sample_rate, filtered_signal)
    # plt.title("Segnale di Test Filtrato con FIR")
    # plt.xlabel("Tempo (s)")
    # plt.ylabel("Ampiezza")
    # plt.grid(True)
    # plt.show()

