# components/op_amp.py

import numpy as np
from utils.helpers import numerical_jacobian

class OpAmp:
    def __init__(self, gain=1e5, input_resistance=1e6, output_resistance=50.0, saturation_voltage=15.0):
        """
        Modello semplificato di un Amplificatore Operazionale (Op-Amp) ideale con alcune non idealità.
        Implementato come una sorgente di tensione controllata in tensione (VCVS) con
        resistenza di ingresso e uscita e saturazione.

        Args:
            gain (float): Guadagno a ciclo aperto dell'Op-Amp (es. 10^5 per un 741).
            input_resistance (float): Resistenza di ingresso differenziale (Ohm).
            output_resistance (float): Resistenza di uscita (Ohm).
            saturation_voltage (float): Tensione di saturazione dell'uscita (Volt).
                                      L'uscita non può superare +/- questo valore.
        """
        if gain <= 0:
            raise ValueError("Il guadagno dell'Op-Amp deve essere positivo.")
        if input_resistance <= 0:
            raise ValueError("La resistenza di ingresso dell'Op-Amp deve essere positiva.")
        if output_resistance < 0:
            raise ValueError("La resistenza di uscita dell'Op-Amp deve essere non negativa.")
        if saturation_voltage <= 0:
            raise ValueError("La tensione di saturazione deve essere positiva.")

        self.gain = float(gain)
        self.Rin = float(input_resistance)
        self.Rout = float(output_resistance)
        self.Vsat = float(saturation_voltage)

        # Stati interni per la Jacobiana se necessari, anche se il modello VCVS è più semplice
        self._input_diff_voltage = 0.0 # V+ - V-
        self._output_voltage = 0.0

    def calculate_output_voltage(self, v_plus, v_minus):
        """
        Calcola la tensione di uscita dell'Op-Amp.
        Vout = Gain * (V+ - V-)
        L'uscita è limitata dalle tensioni di saturazione.
        """
        diff_voltage = v_plus - v_minus
        
        # Guadagno a ciclo aperto
        output_ideal = self.gain * diff_voltage
        
        # Saturazione dell'uscita
        output_saturated = np.clip(output_ideal, -self.Vsat, self.Vsat)
        
        self._input_diff_voltage = diff_voltage # Aggiorna lo stato per conduttanze
        self._output_voltage = output_saturated
        
        return output_saturated

    def calculate_input_current(self, v_plus, v_minus):
        """
        Calcola la corrente che fluisce nei pin di ingresso.
        Per un modello con Rin, la corrente di ingresso è V_diff / Rin.
        """
        diff_voltage = v_plus - v_minus
        # La corrente che entra in V+ e V- è differenziale (se il modello è differenziale)
        # Per un modello ideale, le correnti di ingresso sono zero.
        # Per un modello con Rin: I_in = (V+ - V-) / Rin. Per MNA, è gestito come un resistore tra V+ e V-.
        return diff_voltage / self.Rin # Questa è la corrente che scorre tra i due ingressi.

    # --- Metodi per la Jacobiana ---
    # Questi calcolano le derivate parziali necessarie per la matrice Jacobiana.
    # Per un Op-Amp modellato come VCVS, le conduttanze sono legate al suo comportamento.

    def get_output_resistance(self):
        """Restituisce la resistenza di uscita per il modello VCVS."""
        return self.Rout

    def get_input_resistance(self):
        """Restituisce la resistenza di ingresso differenziale."""
        return self.Rin
        
    def calculate_transconductance_v_plus(self, v_plus, v_minus):
        """
        Calcola la transconduttanza d(Vout)/d(V_plus).
        Questa è semplicemente il guadagno, tenendo conto della saturazione.
        """
        # La derivata di clip è 0 quando è in saturazione, altrimenti è il guadagno.
        if abs(self.gain * (v_plus - v_minus)) >= self.Vsat:
            return 0.0 # In saturazione, la tensione di uscita non cambia con l'input
        return self.gain

    def calculate_transconductance_v_minus(self, v_plus, v_minus):
        """
        Calcola la transconduttanza d(Vout)/d(V_minus).
        Questa è semplicemente -guadagno, tenendo conto della saturazione.
        """
        if abs(self.gain * (v_plus - v_minus)) >= self.Vsat:
            return 0.0
        return -self.gain

    def __str__(self):
        return f"OpAmp(Gain={self.gain:.1e}, Rin={self.Rin:.1e} Ohm, Rout={self.Rout:.1f} Ohm, Vsat={self.Vsat:.1f}V)"
