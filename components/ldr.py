import numpy as np
from utils.constants import DEFAULT_SAMPLE_RATE

class LDR:
    """
    Modello fisico per una fotoresistenza (LDR - Light Dependent Resistor).
    La sua resistenza dipende dall'intensità luminosa incidente e ha una dinamica di risposta.
    """
    def __init__(self, resistance_dark=1e6, resistance_light=1e3,
                 light_reference=1.0, gamma_factor=0.7,
                 time_constant_rise_s=0.01, time_constant_fall_s=0.05,
                 sample_rate=DEFAULT_SAMPLE_RATE):
        """
        Inizializza la LDR.

        Args:
            resistance_dark (float): Resistenza della LDR al buio (Ohm).
            resistance_light (float): Resistenza della LDR sotto una luce di riferimento (Ohm).
            light_reference (float): Valore di luce corrispondente a resistance_light.
                                     Questo è un punto di calibrazione.
            gamma_factor (float): Esponente che descrive la sensibilità della LDR alla luce.
                                  (R = R_dark * (Light_ref / Light_inc)^gamma)
            time_constant_rise_s (float): Costante di tempo di salita (attacco) in secondi.
            time_constant_fall_s (float): Costante di tempo di caduta (rilascio) in secondi.
            sample_rate (float): Frequenza di campionamento in Hz.
        """
        if resistance_dark <= 0 or resistance_light <= 0 or light_reference <= 0:
            raise ValueError("Le resistenze e la luce di riferimento devono essere positive.")
        if gamma_factor <= 0:
            raise ValueError("Il fattore gamma deve essere positivo.")
        if time_constant_rise_s <= 0 or time_constant_fall_s <= 0:
            raise ValueError("Le costanti di tempo devono essere positive.")
        if sample_rate <= 0:
            raise ValueError("La frequenza di campionamento deve essere positiva.")

        self.R_dark = float(resistance_dark)
        self.R_light = float(resistance_light)
        self.L_ref = float(light_reference)
        self.gamma = float(gamma_factor)
        self.Ts = 1.0 / float(sample_rate)

        # Calcola i coefficienti dei filtri passa-basso per le costanti di tempo
        # Implementazione di un filtro passa-basso di primo ordine (alfa per smoothing)
        self.alpha_rise = 1.0 - np.exp(-self.Ts / time_constant_rise_s)
        self.alpha_fall = 1.0 - np.exp(-self.Ts / time_constant_fall_s)

        self.smoothed_light_level = 0.0 # Livello di luce effettivo percepito dalla LDR
        self.current_resistance = self.R_dark # Resistenza attuale della LDR

        print(f"LDR creata (R_dark={self.R_dark} Ohm, gamma={self.gamma}, Rise={time_constant_rise_s*1e3:.0f}ms, Fall={time_constant_fall_s*1e3:.0f}ms).")

    def get_resistance(self, incident_light_level):
        """
        Calcola la resistenza attuale della LDR basata sul livello di luce incidente
        e l'inerzia interna del componente.

        Args:
            incident_light_level (float): Livello di luce (ad esempio, dall'output di un LED).

        Returns:
            float: La resistenza attuale della LDR in Ohm.
        """
        # Aggiorna il livello di luce percepito con un filtro dinamico (attacco/rilascio)
        if incident_light_level > self.smoothed_light_level:
            self.smoothed_light_level += self.alpha_rise * (incident_light_level - self.smoothed_light_level)
        else:
            self.smoothed_light_level += self.alpha_fall * (incident_light_level - self.smoothed_light_level)

        # Evita divisione per zero e log di zero se light_level è molto basso o zero
        # Assicurati che smoothed_light_level sia sempre > 0 per la formula logaritmica
        effective_light = max(self.smoothed_light_level, 1e-10) # Previene lo zero

        # Formula per la resistenza, spesso approssimata da R = R_dark * (L_ref / L_incident)^gamma
        # Possiamo linearizzare per il punto di calibrazione:
        # ln(R) = ln(R_dark) - gamma * (ln(L_inc) - ln(L_ref))
        # Oppure R = R_dark * (L_ref / L_inc)^gamma
        # Un modo comune è interpolare tra R_dark e R_light in scala logaritmica

        # Un modello più robusto che usa il gamma:
        if effective_light <= 0: # Caso di buio completo o assenza di luce
            self.current_resistance = self.R_dark
        else:
            # Assicurati che l'esponente sia calcolabile e non produca valori enormi o infinitesimi
            ratio = self.L_ref / effective_light
            # Limita il rapporto per evitare overflow/underflow con l'esponente
            ratio = np.clip(ratio, 1e-10, 1e10) # Evita valori estremi per l'esponenziazione
            self.current_resistance = self.R_dark * (ratio ** self.gamma)

        # Assicurati che la resistenza sia all'interno di un range sensato
        self.current_resistance = np.clip(self.current_resistance, self.R_light, self.R_dark)

        return self.current_resistance

    def __str__(self):
        return (f"LDR(R_dark={self.R_dark:.1e} Ohm, R_light={self.R_light:.1e} Ohm, "
                f"gamma={self.gamma}, Rise={self.alpha_rise:.3f}, Fall={self.alpha_fall:.3f})")

# Esempio di utilizzo (simulazione della risposta a un impulso di luce)
if __name__ == "__main__":
    fs = DEFAULT_SAMPLE_RATE
    ldr1 = LDR(sample_rate=fs, time_constant_rise_s=0.01, time_constant_fall_s=0.1)
    print(ldr1)

    time_vector = np.arange(0, 1.0, ldr1.Ts) # 1 secondo di simulazione

    # Simuliamo un impulso di luce
    light_input = np.zeros_like(time_vector)
    light_input[int(0.1 * fs):int(0.3 * fs)] = 100.0 # Luce accesa
    light_input[int(0.5 * fs):int(0.7 * fs)] = 50.0 # Luce accesa a metà intensità

    resistance_output = []
    for light_level in light_input:
        resistance_output.append(ldr1.get_resistance(light_level))

    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.plot(time_vector, light_input, label='Luce Incidente')
    plt.title('Livello di Luce Incidente')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Luminosità')
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.plot(time_vector, resistance_output, label='Resistenza LDR')
    plt.title('Resistenza LDR nel Tempo')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Resistenza (Ohm)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
