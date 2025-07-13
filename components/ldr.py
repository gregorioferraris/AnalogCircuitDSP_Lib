# components/ldr.py

import numpy as np

class LDR:
    def __init__(self, R_dark=1e6, R_light=100, light_level_ref=100, alpha=-0.7,
                 sample_rate=48000, time_constant_rise_s=0.01, time_constant_fall_s=0.1):
        """
        Simula una PhotoResistor (Light Dependent Resistor - LDR).
        Args:
            R_dark (float): Resistenza in assenza di luce (Ohm).
            R_light (float): Resistenza a un certo livello di luce di riferimento (Ohm).
            light_level_ref (float): Livello di luce di riferimento per R_light.
            alpha (float): Esponente della curva di risposta R = R0 * (Light/Light0)^alpha.
            sample_rate (float): Frequenza di campionamento per la simulazione dinamica.
            time_constant_rise_s (float): Tempo di salita (tau) in secondi quando la luce aumenta.
            time_constant_fall_s (float): Tempo di caduta (tau) in secondi quando la luce diminuisce.
        """
        self.R_dark = float(R_dark)
        self.R_light = float(R_light)
        self.light_level_ref = float(light_level_ref)
        self.alpha = float(alpha) # Generalmente negativo

        self.sample_rate = float(sample_rate)
        self.Ts = 1.0 / self.sample_rate

        self.time_constant_rise_s = float(time_constant_rise_s)
        self.time_constant_fall_s = float(time_constant_fall_s)

        # Stato interno per la simulazione dinamica della resistenza
        self.current_resistance = R_dark # Inizia al buio
        self.prev_light_level = 0.0

        # Calcola R0 per la formula R = R0 * Light^alpha
        # R_light = R0 * (light_level_ref)^alpha => R0 = R_light / (light_level_ref)^alpha
        self.R0 = self.R_light / (self.light_level_ref ** self.alpha)


    def _calculate_steady_state_resistance(self, light_level):
        """Calcola la resistenza dell'LDR in stato stazionario per un dato livello di luce."""
        if light_level <= 0:
            return self.R_dark
        # Utilizza la formula R = R0 * Light^alpha
        # R = R_light * (light_level / light_level_ref)^alpha
        # Per evitare overflow con np.power se alpha è negativo e light_level è molto piccolo
        light_level_clamped = max(1e-9, light_level) # Prevenire divisione per zero o log(0)
        return self.R0 * (light_level_clamped ** self.alpha)


    def get_resistance(self, light_level):
        """
        Restituisce la resistenza attuale dell'LDR, aggiornata dinamicamente
        in base al livello di luce. Implementa una risposta dinamica.
        """
        target_resistance = self._calculate_steady_state_resistance(light_level)

        # Applica costante di tempo per salita o discesa
        if light_level > self.prev_light_level: # Luce in aumento, resistenza in diminuzione
            tau = self.time_constant_rise_s
        else: # Luce in diminuzione, resistenza in aumento
            tau = self.time_constant_fall_s

        # Filtro passa-basso di primo ordine per modellare la risposta dinamica
        # new_value = old_value + (target_value - old_value) * (1 - exp(-Ts / tau))
        alpha_filter = 1.0 - np.exp(-self.Ts / tau)
        self.current_resistance = self.current_resistance + \
                                  (target_resistance - self.current_resistance) * alpha_filter

        self.prev_light_level = light_level # Aggiorna il livello di luce precedente

        # Assicura che la resistenza non sia mai negativa o zero
        return max(self.current_resistance, 1.0) # Almeno 1 Ohm per evitare divisioni per zero


    def __str__(self):
        return f"LDR(R_dark={self.R_dark:.1e}, R_light={self.R_light:.1f} @ {self.light_level_ref} units)"
