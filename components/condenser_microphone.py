# components/condenser_microphone.py
import numpy as np
from components.component import Component

class CondenserMicrophone(Component):
    def __init__(self, name: str, output_node: str, ground_node: str,
                 diaphragm_area: float = 1e-4, diaphragm_mass: float = 1e-5, diaphragm_compliance: float = 1e-7,
                 polarization_voltage: float = 48.0, capsule_capacitance: float = 50e-12, preamp_input_impedance: float = 1e9):
        """
        Inizializza un modello semplificato di microfono a condensatore.
        Converte la pressione acustica in un segnale elettrico.
        Args:
            name (str): Nome univoco dell'istanza (es. "CMIC1").
            output_node (str): Nodo di uscita del segnale elettrico.
            ground_node (str): Nodo di riferimento (ground).
            diaphragm_area (float): Area effettiva del diaframma (m^2).
            diaphragm_mass (float): Massa del diaframma (kg).
            diaphragm_compliance (float): Cedevolezza meccanica del diaframma (m/N).
            polarization_voltage (float): Tensione di polarizzazione della capsula (Volt).
            capsule_capacitance (float): Capacità nominale della capsula (Farad).
            preamp_input_impedance (float): Impedenza di ingresso del preamplificatore interno (Ohm).
        """
        super().__init__(name, output_node, ground_node)
        self.pin_names = ('output', 'ground')

        self.diaphragm_area = diaphragm_area
        self.diaphragm_mass = diaphragm_mass
        self.diaphragm_compliance = diaphragm_compliance
        self.polarization_voltage = polarization_voltage
        self.capsule_capacitance = capsule_capacitance
        self.preamp_input_impedance = preamp_input_impedance # Importante per la risposta in bassa frequenza

        # Parametri meccanici equivalenti (analogia elettrica)
        # Massa -> Induttanza, Cedevolezza -> Capacità, Resistenza meccanica (damping) -> Resistenza
        # Per semplicità, non modelliamo esplicitamente la resistenza meccanica qui, ma è presente.
        self.L_mech_eq = self.diaphragm_mass
        self.C_mech_eq = self.diaphragm_compliance

    def calculate_output_voltage(self, acoustic_pressure: float) -> float:
        """
        Calcola la tensione di uscita del microfono data la pressione acustica.
        Questo è un modello molto semplificato.
        Args:
            acoustic_pressure (float): Pressione sonora in Pascal (Pa).
        Returns:
            float: Tensione di uscita in Volt.
        """
        # Forza sul diaframma: F = P * A
        force = acoustic_pressure * self.diaphragm_area

        # Spostamento del diaframma (semplificato, ignora la dinamica risonante per ora)
        # x = F * C_mech_eq (per un sistema statico)
        # Per un modello dinamico, la velocità del diaframma sarebbe l'output di un sistema massa-molla-smorzatore.
        # Qui, usiamo un guadagno di sensibilità semplificato.
        
        # Sensibilità del microfono a condensatore (V/Pa)
        # Sensibilità = (Polarization_Voltage / Distanza_Diaframma_Piastra_Fissa) * Area_Diaframma * Cedevolezza
        # Qui, usiamo un fattore di sensibilità approssimato.
        
        # Un modello più accurato calcolerebbe la variazione di capacità e la tensione risultante.
        # Variazione di capacità dC = C0 * (dx / d0)
        # Tensione dV = Vbias * (dC / C0) = Vbias * (dx / d0)
        # dx è lo spostamento, d0 è la distanza a riposo.
        
        # Per ora, un modello lineare basato sulla sensibilità:
        # Assumiamo una sensibilità tipica (es. 10 mV/Pa = 0.01 V/Pa)
        # La sensibilità reale dipende da tutti i parametri fisici.
        
        # Un approccio migliore per il pre-calcolo sarebbe una funzione di trasferimento
        # nel dominio della frequenza, simile agli altoparlanti.
        
        # Placeholder per la sensibilità:
        sensitivity_factor = self.polarization_voltage * self.diaphragm_area * self.diaphragm_compliance * 1e3 # Esempio di scala
        
        output_voltage = acoustic_pressure * sensitivity_factor
        
        # Aggiungi qui eventuali limitazioni o rumore se necessario
        return output_voltage

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Per i microfoni, get_stamps non contribuisce direttamente alla matrice MNA
        nel modo tradizionale, poiché sono generatori di segnale basati su un input acustico.
        La loro uscita è una sorgente di tensione (o corrente) controllata esternamente.
        """
        return np.zeros((num_total_equations, num_total_equations)), np.zeros(num_total_equations)

    def get_frequency_response(self, frequencies: np.ndarray) -> np.ndarray:
        """
        Calcola la risposta in frequenza del microfono a condensatore.
        Questo è un modello semplificato che include la risonanza meccanica del diaframma
        e l'effetto dell'impedenza di ingresso del preamplificatore.
        Args:
            frequencies (np.ndarray): Frequenze in Hz.
        Returns:
            np.ndarray: Risposta in frequenza complessa (V/Pa).
        """
        omega = 2 * np.pi * frequencies
        
        # Risonanza meccanica del diaframma (massa-molla-smorzatore)
        # Z_mech = R_mech + jwM_mech + 1/(jwC_mech)
        # Per ora, assumiamo R_mech sia piccolo o incluso in un fattore Q.
        # Frequenza di risonanza F_res = 1 / (2*pi*sqrt(L_mech_eq * C_mech_eq))
        
        # La funzione di trasferimento dalla forza alla velocità è 1/Z_mech
        # La velocità del diaframma è proporzionale alla tensione di uscita.
        
        # Un modello semplificato della risposta in frequenza:
        # Un filtro passa-alto dovuto all'impedenza di ingresso del preamplificatore (se AC coupled)
        # e una risonanza di secondo ordine per il diaframma.
        
        # Impedenza meccanica del diaframma (analogia elettrica)
        Z_mech_diaphragm = 1j * omega * self.L_mech_eq + 1.0 / (1j * omega * self.C_mech_eq)
        # Aggiungere un termine di resistenza meccanica (damping) se necessario
        
        # La sensibilità è proporzionale a 1 / (Z_mech_diaphragm) per la velocità
        # e poi alla variazione di capacità per la tensione.
        
        # Per ora, un filtro passa-alto RC per l'accoppiamento con il preamplificatore
        # e una risonanza semplificata.
        
        # Filtro passa-alto RC dovuto alla capsula e all'impedenza di ingresso del preamp
        # F_cut = 1 / (2 * pi * R_preamp * C_capsule)
        Z_capsule = 1.0 / (1j * omega * self.capsule_capacitance)
        Z_preamp = self.preamp_input_impedance # Resistenza di ingresso del preamp
        
        # La tensione di uscita è proporzionale alla corrente che fluisce attraverso Z_preamp
        # I_capsule = V_polarization * jw * C_capsule * (dx / d0)
        # V_out = I_capsule * Z_preamp
        
        # Questo è un modello molto semplificato per la risposta in frequenza.
        # Per una risposta più realistica, si dovrebbe modellare l'accoppiamento acustico
        # e la dinamica completa del diaframma.
        
        # Per ora, restituiamo un guadagno costante con una risonanza semplificata.
        # Questo è un placeholder per la funzione di trasferimento reale.
        
        # Esempio: un picco di risonanza e un roll-off in bassa frequenza
        # Freq. di risonanza meccanica: fres = 1 / (2 * pi * sqrt(L_mech_eq * C_mech_eq))
        # Q-factor del risonatore meccanico
        
        # Per una risposta flat, restituisci un valore costante.
        # Per modellare la risposta reale, dovresti derivare la funzione di trasferimento
        # dal modello fisico completo del microfono.
        
        # Esempio di risposta semplificata (flat con roll-off in bassa frequenza)
        # Questo è un filtro passa-alto di primo ordine.
        # La frequenza di taglio dipende da C_capsule e R_preamp.
        Rc = self.preamp_input_impedance
        Cc = self.capsule_capacitance
        
        # H(jw) = (jwRcCc) / (1 + jwRcCc)
        # Questo è solo il filtro di accoppiamento, non la risposta acustica.
        # La sensibilità acustica dovrebbe essere un fattore moltiplicativo.
        
        # Per un modello più realistico, si dovrebbe considerare il guadagno di sensibilità
        # e la risposta meccanica.
        
        # Placeholder: guadagno costante per la sensibilità, con un roll-off in bassa frequenza.
        # Esempio: 0.01 V/Pa (10 mV/Pa)
        sensitivity_gain = 0.01 # V/Pa
        
        # Filtro passa-alto RC per l'accoppiamento (per modellare la risposta in bassa frequenza)
        cutoff_freq_rc = 1.0 / (2 * np.pi * Rc * Cc)
        
        response = sensitivity_gain * (1j * omega * Rc * Cc) / (1 + 1j * omega * Rc * Cc)
        
        # Aggiungi qui la risonanza meccanica se desideri
        # Esempio di risonanza (filtro passa-banda di secondo ordine):
        # f_res = 1.0 / (2 * np.pi * np.sqrt(self.L_mech_eq * self.C_mech_eq))
        # Q = ... (dipende dallo smorzamento)
        # H_res = (1j * omega * f_res / Q) / ( (1j * omega)**2 + (1j * omega * f_res / Q) + f_res**2 )
        # response *= H_res # Moltiplica per la risposta risonante
        
        return response

