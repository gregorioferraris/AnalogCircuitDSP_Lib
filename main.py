# main.py - Script di test ed esempio per la libreria di simulazione circuitale

import numpy as np
import matplotlib.pyplot as plt

# Importa tutti i componenti implementati finora
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor
from components.diode import Diode
from components.led import LED
from components.ldr import LDR
from components.jfet import JFET
from components.bjt import BJT
from components.schottky_diode import SchottkyDiode
from components.zener_diode import ZenerDiode
from components.mosfet import MOSFET
from components.triode import Triode
from components.pentode import Pentode
from components.rectifier_tube import RectifierTube

from utils.constants import DEFAULT_SAMPLE_RATE # Importa la frequenza di campionamento predefinita


print("--- Inizio Test dei Componenti Discreti ---")
print(f"Frequenza di Campionamento Predefinita: {DEFAULT_SAMPLE_RATE} Hz")

# --- 1. Test Resistore ---
print("\n--- Testing Resistor ---")
r1 = Resistor(1000) # 1 kOhm
voltage_test_r = 9.0 # 9 Volt
current_r = r1.calculate_current(voltage_test_r)
print(f"Su un {r1}, con {voltage_test_r}V, la corrente è {current_r:.4f} A")


# --- 2. Test Condensatore (comportamento dinamico semplificato) ---
print("\n--- Testing Capacitor (simplified dynamic response) ---")
fs = DEFAULT_SAMPLE_RATE
c1 = Capacitor(1e-6, sample_rate=fs) # 1 uF
time_vector_c = np.arange(0, 0.05, c1.Ts) # 50 ms di simulazione
freq_hz_c = 100 # Frequenza della sinusoide
input_voltage_c = 5.0 * np.sin(2 * np.pi * freq_hz_c * time_vector_c) # 5Vpk sinusoide

output_current_c = []
for v_now_c in input_voltage_c:
    i_now_c = c1.calculate_current(v_now_c)
    output_current_c.append(i_now_c)
    c1.update_state(v_now_c, i_now_c) # Aggiorna lo stato del condensatore

plt.figure(figsize=(10, 6))
plt.plot(time_vector_c, input_voltage_c, label='Tensione Ingresso (V)')
plt.plot(time_vector_c, output_current_c, label='Corrente Uscita (A)')
plt.title(f'Risposta Dinamica del Condensatore ({c1.C*1e6:.0f} uF)')
plt.xlabel('Tempo (s)')
plt.ylabel('Valore')
plt.legend()
plt.grid(True)
plt.show()


# --- 3. Test Induttore (comportamento dinamico semplificato) ---
print("\n--- Testing Inductor (simplified dynamic response) ---")
l1 = Inductor(10e-3, sample_rate=fs) # 10 mH
time_vector_l = np.arange(0, 0.05, l1.Ts) # 50 ms di simulazione
freq_hz_l = 100 # Frequenza della sinusoide
input_voltage_l = 5.0 * np.sin(2 * np.pi * freq_hz_l * time_vector_l) # Stessa sinusoide

output_current_l = []
for v_now_l in input_voltage_l:
    i_now_l = l1.calculate_current(v_now_l)
    output_current_l.append(i_now_l)
    l1.update_state(v_now_l, i_now_l) # Aggiorna lo stato dell'induttore

plt.figure(figsize=(10, 6))
plt.plot(time_vector_l, input_voltage_l, label='Tensione Ingresso (V)')
plt.plot(time_vector_l, output_current_l, label='Corrente Uscita (A)')
plt.title(f'Risposta Dinamica dell\'Induttore ({l1.L*1e3:.0f} mH)')
plt.xlabel('Tempo (s)')
plt.ylabel('Valore')
plt.legend()
plt.grid(True)
plt.show()


# --- 4. Test Diodo Standard ---
print("\n--- Testing Standard Diode ---")
d1 = Diode()
print(d1)
voltages_d = np.linspace(-1, 1, 200) # Da -1V a 1V
currents_d = np.array([d1.calculate_current(v) for v in voltages_d])
plt.figure(figsize=(8, 5))
plt.plot(voltages_d, currents_d * 1e3) # Corrente in mA
plt.title('Caratteristica I-V del Diodo Standard')
plt.xlabel('Tensione (V)')
plt.ylabel('Corrente (mA)')
plt.grid(True)
plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
plt.show()


# --- 5. Test LED ---
print("\n--- Testing LED ---")
led1 = LED(luminous_efficiency=100)
print(led1)
voltages_led = np.linspace(-1, 2, 200) # Da -1V a 2V
currents_led = np.array([led1.calculate_current(v) for v in voltages_led])
luminous_outputs_led = led1.get_luminous_output(currents_led)
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(voltages_led, currents_led * 1e3) # Corrente in mA
plt.title('Caratteristica I-V del LED')
plt.xlabel('Tensione (V)')
plt.ylabel('Corrente (mA)')
plt.grid(True)
plt.subplot(2, 1, 2)
plt.plot(currents_led, luminous_outputs_led)
plt.title('Output Luminoso del LED vs Corrente')
plt.xlabel('Corrente (A)')
plt.ylabel('Luminosità Relativa')
plt.grid(True)
plt.tight_layout()
plt.show()


# --- 6. Test LDR ---
print("\n--- Testing LDR ---")
ldr1 = LDR(sample_rate=fs, time_constant_rise_s=0.01, time_constant_fall_s=0.1)
print(ldr1)
time_vector_ldr = np.arange(0, 1.0, ldr1.Ts) # 1 secondo di simulazione
light_input_ldr = np.zeros_like(time_vector_ldr)
light_input_ldr[int(0.1 * fs):int(0.3 * fs)] = 100.0 # Luce accesa
light_input_ldr[int(0.5 * fs):int(0.7 * fs)] = 50.0 # Luce accesa a metà intensità

resistance_output_ldr = []
for light_level in light_input_ldr:
    resistance_output_ldr.append(ldr1.get_resistance(light_level))

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time_vector_ldr, light_input_ldr, label='Luce Incidente')
plt.title('Livello di Luce Incidente (LDR)')
plt.xlabel('Tempo (s)')
plt.ylabel('Luminosità')
plt.grid(True)
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(time_vector_ldr, resistance_output_ldr, label='Resistenza LDR')
plt.title('Resistenza LDR nel Tempo')
plt.xlabel('Tempo (s)')
plt.ylabel('Resistenza (Ohm)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# --- 7. Test JFET ---
print("\n--- Testing JFET ---")
jfet1 = JFET(Idss=0.01, Vp=-2.0, lambda_val=0.05) # JFET N-channel
print(jfet1)
vds_values_jfet = np.linspace(0, 5, 100) # Vds da 0 a 5V
vgs_steps_jfet = np.linspace(jfet1.Vp, 0, 5) # Diversi valori di Vgs
plt.figure(figsize=(10, 6))
for vgs in vgs_steps_jfet:
    id_values_jfet = jfet1.calculate_drain_current(vgs, vds_values_jfet)
    plt.plot(vds_values_jfet, id_values_jfet * 1e3, label=f'Vgs={vgs:.1f}V') # Corrente in mA
plt.title('Caratteristiche di Uscita del JFET (Id vs Vds)')
plt.xlabel('Vds (V)')
plt.ylabel('Id (mA)')
plt.legend()
plt.grid(True)
plt.show()


# --- 8. Test BJT ---
print("\n--- Testing BJT ---")
bjt1 = BJT(Is=1e-14, Bf=200, Vaf=80.0) # BJT NPN
print(bjt1)
vce_values_bjt = np.linspace(0, 10, 100) # Vce da 0 a 10V
vbe_steps_bjt = np.linspace(0.6, 0.75, 5) # Vbe da 0.6V a 0.75V
plt.figure(figsize=(10, 6))
for vbe in vbe_steps_bjt:
    ic_values_bjt = bjt1.calculate_collector_current(vbe, vce_values_bjt)
    plt.plot(vce_values_bjt, ic_values_bjt * 1e3, label=f'Vbe={vbe:.2f}V') # Corrente in mA
plt.title('Caratteristiche di Uscita del BJT (Ic vs Vce)')
plt.xlabel('Vce (V)')
plt.ylabel('Ic (mA)')
plt.legend()
plt.grid(True)
plt.show()


# --- 9. Test Diodo Zener ---
print("\n--- Testing Zener Diode ---")
zener1 = ZenerDiode(Vz=5.6, Iz_test=1e-3, Rz=10.0) # Zener 5.6V
print(zener1)
voltages_z = np.linspace(-10, 1, 200) # Da -10V a 1V
currents_z = np.array([zener1.calculate_current(v) for v in voltages_z])
plt.figure(figsize=(10, 6))
plt.plot(voltages_z, currents_z * 1e3) # Corrente in mA
plt.title(f'Caratteristica I-V del Diodo Zener ({zener1.Vz}V)')
plt.xlabel('Tensione (V)')
plt.ylabel('Corrente (mA)')
plt.grid(True)
plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
plt.annotate(f'Vz = {-zener1.Vz}V', xy=(-zener1.Vz, 0), xytext=(-zener1.Vz - 2, 5),
             arrowprops=dict(facecolor='black', shrink=0.05),
             horizontalalignment='right')
plt.show()


# --- 10. Test Diodo Schottky ---
print("\n--- Testing Schottky Diode ---")
schottky1 = SchottkyDiode()
std_diode_for_comp = Diode() # Diodo standard per confronto
print(schottky1)
voltages_s = np.linspace(-0.5, 1.0, 200) # Da -0.5V a 1.0V
schottky_currents_s = np.array([schottky1.calculate_current(v) for v in voltages_s])
std_diode_currents_s = np.array([std_diode_for_comp.calculate_current(v) for v in voltages_s])
plt.figure(figsize=(10, 6))
plt.plot(voltages_s, schottky_currents_s * 1e3, label='Diodo Schottky (mA)')
plt.plot(voltages_s, std_diode_currents_s * 1e3, label='Diodo Standard (mA)', linestyle='--')
plt.title('Caratteristica I-V: Schottky vs Diodo Standard')
plt.xlabel('Tensione (V)')
plt.ylabel('Corrente (mA)')
plt.legend()
plt.grid(True)
plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
plt.show()


# --- 11. Test MOSFET ---
print("\n--- Testing MOSFET ---")
mosfet1 = MOSFET(Vt=1.5, Kn=1.0e-3, lambda_val=0.02) # MOSFET N-Enh
print(mosfet1)
vds_values_mosfet = np.linspace(0, 10, 100) # Vds da 0 a 10V
vgs_steps_mosfet = np.linspace(mosfet1.Vt, mosfet1.Vt + 5, 5) # Diversi Vgs sopra la soglia
plt.figure(figsize=(10, 6))
for vgs in vgs_steps_mosfet:
    id_values_mosfet = mosfet1.calculate_drain_current(vgs, vds_values_mosfet)
    plt.plot(vds_values_mosfet, id_values_mosfet * 1e3, label=f'Vgs={vgs:.1f}V') # Corrente in mA
plt.title('Caratteristiche di Uscita del MOSFET (Id vs Vds)')
plt.xlabel('Vds (V)')
plt.ylabel('Id (mA)')
plt.legend()
plt.grid(True)
plt.show()


# --- 12. Test Triodo ---
print("\n--- Testing Triode ---")
triode1 = Triode(mu=100.0, Kp=600.0, X=1.5, Kg1=2.5) # 12AX7 style
print(triode1)
va_values_triode = np.linspace(0, 400, 100) # Va da 0 a 400V
vg_steps_triode = np.array([0, -0.5, -1, -2, -4, -8]) # Diversi Vg
plt.figure(figsize=(10, 6))
for vg in vg_steps_triode:
    ia_values_triode = triode1.calculate_anode_current(vg, va_values_triode)
    plt.plot(va_values_triode, ia_values_triode * 1e3, label=f'Vg={vg:.1f}V') # Corrente in mA
plt.title('Caratteristiche di Uscita del Triodo (Ia vs Va)')
plt.xlabel('Va (V)')
plt.ylabel('Ia (mA)')
plt.legend()
plt.grid(True)
plt.show()


# --- 13. Test Pentodo ---
print("\n--- Testing Pentode ---")
pentode1 = Pentode(mu=10.0, Kp=300.0, X=1.5, Kg1=5.0, Kg2=10.0) # EL34/6L6 style
print(pentode1)
va_values_pentode = np.linspace(0, 400, 100) # Va da 0 a 400V
vg1_steps_pentode = np.array([-5, -10, -15, -20]) # Diversi Vg1
vg2_fixed_pentode = 250.0 # Vg2 fissa
plt.figure(figsize=(10, 6))
for vg1 in vg1_steps_pentode:
    ia_values_pentode = pentode1.calculate_anode_current(vg1, vg2_fixed_pentode, va_values_pentode)
    plt.plot(va_values_pentode, ia_values_pentode * 1e3, label=f'Vg1={vg1:.1f}V') # Corrente in mA
plt.title(f'Caratteristiche di Uscita del Pentodo (Ia vs Va) con Vg2={vg2_fixed_pentode}V')
plt.xlabel('Va (V)')
plt.ylabel('Ia (mA)')
plt.legend()
plt.grid(True)
plt.show()


# --- 14. Test Diodo Raddrizzatore a Valvola ---
print("\n--- Testing Rectifier Tube ---")
rectifier1 = RectifierTube(forward_voltage_drop=15.0, series_resistance=100.0) # Tipico 5AR4
print(rectifier1)
voltages_rect = np.linspace(-50, 50, 200) # Da -50V a 50V
currents_rect = np.array([rectifier1.calculate_current(v) for v in voltages_rect])
plt.figure(figsize=(10, 6))
plt.plot(voltages_rect, currents_rect * 1e3) # Corrente in mA
plt.title('Caratteristica I-V del Diodo Raddrizzatore a Valvola')
plt.xlabel('Tensione (V)')
plt.ylabel('Corrente (mA)')
plt.grid(True)
plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
plt.show()


print("\n--- Tutti i test dei componenti discreti completati. ---")
print("Questi test mostrano le curve caratteristiche individuali. Il prossimo passo è integrarli in un solutore di circuito per simulare intere reti.")
