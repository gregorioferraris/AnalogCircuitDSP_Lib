import numpy as np

# Costanti fisiche
PI = np.pi               # Pi greco
K = 1.380649e-23         # Costante di Boltzmann (J/K)
Q = 1.602176634e-19      # Carica elementare (Coulomb)
T_ROOM = 298.15          # Temperatura ambiente in Kelvin (25 °C)
VT_ROOM = K * T_ROOM / Q # Tensione termica a temperatura ambiente (circa 25.7mV)

# Parametri di simulazione globali (possono essere sovrascritti dal solutore)
DEFAULT_SAMPLE_RATE = 44100.0 # Hz (Frequenza di campionamento standard)
