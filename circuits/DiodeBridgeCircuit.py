# circuits/DiodeBridgeCompressorSidechain.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from circuits.DiodeBridgeCircuit import DiodeBridgeCircuit # Importiamo il ponte di diodi

import numpy as np

class DiodeBridgeCompressorSidechain(Circuit):
    """
    Circuito completo per la side-chain di un compressore basato su ponte di diodi.
    Include:
    1. Un ponte di diodi per rettificare il segnale audio (o copia del segnale).
    2. Un circuito RC per l'envelope detection, che determina i tempi di Attack e Release.
    L'output di questo circuito è la tensione di controllo DC per il VCA (o elemento di guadagno).
    """
    def __init__(self, name="DiodeSidechain",
                 diode_is=1e-14, diode_n=1.0,
                 attack_resistor=10000.0, release_capacitor=1.0e-6, # RC per attacco/rilascio
                 load_resistor=100000.0, # Resistenza di carico per il condensatore di release
                 sample_rate=48000):
        super().__init__(name)
        self.sample_rate = sample_rate

        self.attack_R_val = float(attack_resistor)
        self.release_C_val = float(release_capacitor)
        self.load_R_val = float(load_resistor) # Aggiungo un carico esplicito per la scarica

        # Inizializza il ponte di diodi come sottocircuito
        self.diode_bridge = DiodeBridgeCircuit(name="InternalDiodeBridge",
                                                diode_is=diode_is, diode_n=diode_n,
                                                filter_capacitor_value=self.release_C_val, # Usiamo il C di release qui
                                                load_resistor_value=self.load_R_val,
                                                sample_rate=sample_rate)
        
        # Aggiungi il sottocircuito al circuito principale.
        # Nota: l'MNA dovrà sapere come gestire i nodi dei sottocircuiti.
        # Per semplicità qui, "inglobiamo" il ponte e re-mappiamo i suoi nodi.
        # Un approccio più robusto potrebbe prevedere un meccanismo di aggiunta sottocircuito nel solver.
        
        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()

    def _add_nodes(self):
        """Aggiunge i nodi specifici per la side-chain."""
        self.add_node("Audio_Input") # Ingresso del segnale audio per la side-chain
        self.add_node("Control_Voltage_Output") # Uscita della tensione di controllo DC

        # Nodi interni per il ponte di diodi
        self.add_node("DiodeBridge_AC_Input_P")
        self.add_node("DiodeBridge_AC_Input_N")
        self.add_node("DiodeBridge_DC_Output_P")
        self.add_node("DiodeBridge_DC_Output_N")

        # Nodo tra la resistenza di attack e il condensatore di release
        self.add_node("Envelope_RC_Node")


    def _add_components(self):
        """Aggiunge i componenti: ponte di diodi e filtro RC per l'inviluppo."""
        
        # --- Componenti del Ponte di Diodi (modellati qui direttamente per semplicità di integrazione) ---
        # Avremmo dovuto gestire un sottocircuito, ma per MNA semplifichiamo.
        # Diodi del ponte:
        D1 = self.diode_bridge.diodes[0] # Uso i diodi già configurati nel sottocircuito per i parametri
        D2 = self.diode_bridge.diodes[1]
        D3 = self.diode_bridge.diodes[2]
        D4 = self.diode_bridge.diodes[3]

        self.add_component(D1, "DiodeBridge_AC_Input_P", "DiodeBridge_DC_Output_P")
        self.add_component(D2, "DiodeBridge_AC_Input_N", "DiodeBridge_DC_Output_P")
        self.add_component(D3, "DiodeBridge_DC_Output_N", "DiodeBridge_AC_Input_P")
        self.add_component(D4, "DiodeBridge_DC_Output_N", "DiodeBridge_AC_Input_N")
        
        # --- Circuito RC per l'Envelope Detection (Attack/Release) ---
        # La resistenza di attack è in serie all'uscita del ponte
        self.add_component(Resistor(self.attack_R_val), "DiodeBridge_DC_Output_P", "Envelope_RC_Node")
        
        # Il condensatore di release è in parallelo alla resistenza di carico
        self.add_component(Capacitor(self.release_C_val, sample_rate=self.sample_rate), "Envelope_RC_Node", "GND")
        
        # La resistenza di carico è in parallelo al condensatore di release (importante per il tempo di release)
        self.add_component(Resistor(self.load_R_val), "Envelope_RC_Node", "GND")

    def _connect_nodes(self):
        """Connette i nodi del circuito."""
        # L'ingresso audio va all'ingresso del ponte di diodi
        self.connect_nodes("Audio_Input", "DiodeBridge_AC_Input_P")
        self.connect_nodes("DiodeBridge_AC_Input_N", "GND") # L'altro capo del ponte a GND

        # L'uscita della side-chain è il nodo del condensatore di release
        self.connect_nodes("Control_Voltage_Output", "Envelope_RC_Node")

    def get_audio_input_node(self):
        return self.get_node_id("Audio_Input")

    def get_control_voltage_output_node(self):
        return self.get_node_id("Control_Voltage_Output")

    def get_attack_time(self):
        """
        Calcola il tempo di attack.
        Per un filtro RC, il tempo di carica (attack) è dato da R_attack * C_release.
        """
        return self.attack_R_val * self.release_C_val

    def get_release_time(self):
        """
        Calcola il tempo di release.
        Per un filtro RC, il tempo di scarica (release) è dato da R_load * C_release.
        """
        return self.load_R_val * self.release_C_val

    def set_attack_time(self, attack_ms):
        """Imposta il tempo di attack modificando la resistenza di attack."""
        self.attack_R_val = (attack_ms / 1000.0) / self.release_C_val
        # Aggiorna il componente Resistore nel circuito (necessita un meccanismo nel solver)
        # Per ora, è solo un aggiornamento del valore interno.
        print(f"Tempo di Attack impostato a {attack_ms}ms, R_attack = {self.attack_R_val:.2f} Ohm")

    def set_release_time(self, release_ms):
        """Imposta il tempo di release modificando la resistenza di carico."""
        self.load_R_val = (release_ms / 1000.0) / self.release_C_val
        # Aggiorna il componente Resistore nel circuito (necessita un meccanismo nel solver)
        print(f"Tempo di Release impostato a {release_ms}ms, R_load = {self.load_R_val:.2f} Ohm")
