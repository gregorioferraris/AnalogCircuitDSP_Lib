# main.py - Script di test e esempio per la libreria di simulazione circuitale

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor
from utils.constants import DEFAULT_SAMPLE_RATE # Importiamo la frequenza di campionamento predefinita

import numpy as np
import matplotlib.pyplot as plt # Utile per visualizzare i risultati

print("--- Test dei Componenti Base ---")

# --- Test Resistore ---
print("\nTesting Resistor:")
r1 = Resistor(1000) # 1 kOhm
voltage_test_r = 9.0 # 9 Volt
current_r = r1.calculate_current(voltage_test_r)
print(f"Su un {r1}, con {voltage_test_r}V, la corrente è {current_r:.4f} A")

# --- Test Condensatore (semplificato, richiede un solutore per essere realistico) ---
print("\nTesting Capacitor (simplified, needs solver for full realism):")
fs = DEFAULT_SAMPLE_RATE
c1 = Capacitor(1e-6, sample_rate=fs) # 1 uF capacitor

# Simuliamo una tensione sinusoidale ai capi del condensatore
time_vector = np.arange(0, 0.05, c1.Ts) # 50 ms di simulazione
freq_hz = 100 # Frequenza della sinusoide
input_voltage_c = 5.0 * np.sin(2 * np.pi * freq_hz * time_vector) # 5Vpk sinusoide

output_current_c = []
for v_now_c in input_voltage_c:
    i_now_c = c1.calculate_current(v_now_c)
    output_current_c.append(i_now_c)
    c1.update_state(v_now_c, i_now_c) # Aggiorna lo stato del condensatore

# Plot dei risultati del condensatore
plt.figure(figsize=(10, 6))
plt.plot(time_vector, input_voltage_c, label='Input Voltage (V)')
plt.plot(time_vector, output_current_c, label='Output Current (A)')
plt.title(f'Capacitor ({c1.C*1e6:.0f} uF) Response to Sine Wave')
plt.xlabel('Time (s)')
plt.ylabel('Value')
plt.legend()
plt.grid(True)
plt.show()


# --- Test Induttore (semplificato, richiede un solutore per essere realistico) ---
print("\nTesting Inductor (simplified, needs solver for full realism):")
l1 = Inductor(10e-3, sample_rate=fs) # 10 mH inductor

# Simuliamo una tensione sinusoidale ai capi dell'induttore
input_voltage_l = 5.0 * np.sin(2 * np.pi * freq_hz * time_vector) # Stessa sinusoide

output_current_l = []
for v_now_l in input_voltage_l:
    i_now_l = l1.calculate_current(v_now_l)
    output_current_l.append(i_now_l)
    l1.update_state(v_now_l, i_now_l) # Aggiorna lo stato dell'induttore

# Plot dei risultati dell'induttore
plt.figure(figsize=(10, 6))
plt.plot(time_vector, input_voltage_l, label='Input Voltage (V)')
plt.plot(time_vector, output_current_l, label='Output Current (A)')
plt.title(f'Inductor ({l1.L*1e3:.0f} mH) Response to Sine Wave')
plt.xlabel('Time (s)')
plt.ylabel('Value')
plt.legend()
plt.grid(True)
plt.show()

print("\n--- Test dei componenti base completato ---")
print("I risultati di Capacitore e Induttore sono semplificati, il loro comportamento completo si vedrà con il solutore di circuito.")
