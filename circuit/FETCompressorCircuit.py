# circuits/FETCompressorCircuit.py

from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.fet import FET # Il nostro nuovo componente FET
from components.op_amp import OpAmp # Se usiamo Op-Amp nel percorso audio
from circuits.DiodeBridgeCompressorSidechain import DiodeBridgeCompressorSidechain

import numpy as np

class FETCompressorCircuit(Circuit):
    """
    Circuito di un compressore audio basato su FET.
    Incorpora un FET come elemento di guadagno variabile e una side-chain
    per generare la tensione di controllo.

    Componenti chiave:
    - FET: agisce come resistore variabile (VCR) nel percorso del segnale audio.
    - DiodeBridgeCompressorSidechain: rileva l'inviluppo del segnale per la tensione di controllo.
    - Rete di controllo: collega la tensione di inviluppo al gate del FET.
    - Circuito VCA: in questo caso, il FET stesso inserito in una configurazione di VCA.
    """
    def __init__(self, name="FET_Compressor",
                 # Parametri del FET
                 fet_vt=-2.0, fet_kp=0.005, fet_rds_on_at_zero_vgs=100.0,
                 # Parametri della side-chain (ponte diodi + RC)
                 sidechain_diode_is=1e-14, sidechain_diode_n=1.0,
                 sidechain_attack_ms=10.0, sidechain_release_ms=100.0,
                 # Parametri del circuito VCA (con il FET)
                 vca_input_resistor=10000.0, vca_feedback_resistor=20000.0,
                 sample_rate=48000):
        super().__init__(name)
        self.sample_rate = sample_rate

        self.fet_params = {
            "vt": fet_vt,
            "kp": fet_kp,
            "channel_resistance_at_zero_vgs": fet_rds_on_at_zero_vgs
        }
        self.sidechain_params = {
            "diode_is": sidechain_diode_is,
            "diode_n": sidechain_diode_n,
            "attack_resistor": sidechain_attack_ms, # useremo questi per calcolare R e C
            "release_capacitor": sidechain_release_ms,
            "load_resistor": None, # calcolato da release_ms
            "sample_rate": sample_rate
        }
        self.vca_R_in = float(vca_input_resistor)
        self.vca_R_fb = float(vca_feedback_resistor)

        print(f"Costruendo il circuito: {self.name}")
        self._add_nodes()
        self._add_components()
        self._connect_nodes()
        self._setup_sidechain_times(sidechain_attack_ms, sidechain_release_ms)

    def _add_nodes(self):
        """Aggiunge i nodi per l'intero compressore."""
        self.add_node("Audio_Input")
        self.add_node("Audio_Output")

        # Nodi per la side-chain (saranno connessi al sottocircuito)
        self.add_node("SC_Audio_Input")
        self.add_node("SC_Control_Voltage_Output")

        # Nodi per il FET e il circuito VCA
        self.add_node("FET_Drain")
        self.add_node("FET_Gate")
        self.add_node("FET_Source")

        self.add_node("VCA_OpAmp_Vminus") # Ingresso invertente dell'Op-Amp VCA
        self.add_node("VCA_OpAmp_Vplus")  # Ingresso non invertente dell'Op-Amp VCA
        self.add_node("VCA_OpAmp_Vout")   # Uscita dell'Op-Amp VCA


    def _add_components(self):
        """Aggiunge tutti i componenti: FET, side-chain e VCA."""
        # --- Side-chain: Rilevamento dell'inviluppo ---
        # Creiamo un'istanza della side-chain
        # Nota: i parametri R, C della sidechain saranno impostati da set_attack/release_time
        # Inizializziamo con valori base e poi li settiamo correttamente.
        initial_C = 1e-6 # Valore iniziale per il condensatore di release
        initial_R_attack = 10000.0
        initial_R_load = 100000.0

        self.sidechain = DiodeBridgeCompressorSidechain(
            diode_is=self.sidechain_params["diode_is"],
            diode_n=self.sidechain_params["diode_n"],
            attack_resistor=initial_R_attack,
            release_capacitor=initial_C,
            load_resistor=initial_R_load,
            sample_rate=self.sample_rate
        )
        # Aggiungiamo i componenti della side-chain direttamente (per la gestione MNA semplificata)
        # Questo è un po' ripetitivo, ma necessario se non abbiamo un sistema di sottocircuiti robusto.
        
        # DiodeBridge components:
        # Re-aggiungiamo i diodi della sidechain e i suoi resistori/capacitori.
        # Questo assume che il solutore gestisca i nodi tra i circuiti correttamente.
        # ALTERNATIVA: la sidechain dovrebbe aggiungere i suoi componenti al circuito "padre"
        # Per ora, li aggiungo manualmente qui per chiarezza.

        # I nodi della sidechain vanno mappati:
        # Audio_Input -> SC_Audio_Input
        # Control_Voltage_Output -> SC_Control_Voltage_Output
        # E i nodi interni della sidechain devono essere unici.
        # Questo è il punto debole di non avere un sistema di gestione sottocircuiti robusto.
        # Per aggirarlo, trattiamo la sidechain come un "black box" che genera una tensione,
        # e poi la applichiamo al FET.

        # **Simplification for MNA**: Instead of directly inserting all sidechain components,
        # we will assume the sidechain produces a control voltage at SC_Control_Voltage_Output.
        # The MNA will need to calculate this voltage first.
        # This means the sidechain needs to be solved *before* the main circuit in an iterative manner.
        # For a truly unified MNA, all components must be in one matrix.
        # Let's try to put all components in one circuit:

        # Diodi della Sidechain (mappati a nodi locali del compressore)
        # Questo richiede una profonda ristrutturazione di come DiodeBridgeCircuit è creato.
        # Per mantenere la modularità, assumiamo una "tensione di ingresso" che controlla il FET
        # e calcoliamo quella tensione esternamente al principale solve MNA loop, o facciamo un iterativo.

        # **APPROCCIO PRAGMATICO PER MNA**:
        # La sidechain sarà un circuito separato che viene risolto per produrre il CV.
        # Poi questo CV viene applicato come una sorgente di tensione controllata al Gate del FET.

        # Creiamo il FET
        self.fet = FET("Compressor_FET", **self.fet_params)
        # Setta i nodi del FET
        self.fet.set_nodes(self.get_node_id("FET_Source"),
                           self.get_node_id("FET_Gate"),
                           self.get_node_id("FET_Drain"))
        
        # Non aggiungiamo il FET come "componente" MNA diretto qui,
        # ma i suoi comportamenti I-V e la Jacobiana saranno integrati dal solver.
        self.nonlinear_components.append(self.fet) # Il solutore sa come gestire i componenti non lineari

        # --- VCA Circuit (using an Op-Amp as buffer/inverter and the FET) ---
        # Una configurazione comune è un Op-Amp invertente con il FET nel feedback o in serie.
        # Qui useremo il FET in un divisore di tensione per modulare il segnale.
        # O, più semplicemente, il FET in retroazione per un VCA controllato.

        # Op-Amp per il VCA (guadagno unitario o con guadagno fisso)
        # L'Op-Amp agirà da buffer o inverter, mentre il FET modula il segnale.
        self.vca_op_amp = OpAmp(gain=1e5, input_resistance=1e6, output_resistance=50.0)
        # L'Op-Amp sarà usato in configurazione invertente per il VCA.
        self.add_component(self.vca_op_amp, "VCA_OpAmp_Vplus", "VCA_OpAmp_Vminus", "VCA_OpAmp_Vout")
        
        # Connessione dell'ingresso non invertente dell'Op-Amp a GND
        self.connect_nodes("VCA_OpAmp_Vplus", "GND")

        # Resistenza di ingresso del VCA
        self.add_component(Resistor(self.vca_R_in), "Audio_Input", "VCA_OpAmp_Vminus")

        # Configurazione del FET nel VCA:
        # Il FET è spesso in serie o in parallelo nel percorso di feedback.
        # Qui, mettiamo il FET in serie con la resistenza di feedback per un VCA invertente.
        # La resistenza del FET (Rds) cambia il guadagno: Gain = -(R_fb + Rds_FET) / R_in
        self.add_component(Resistor(self.vca_R_fb), "VCA_OpAmp_Vout", "FET_Drain") # R_fb in serie con FET
        self.connect_nodes("FET_Source", "VCA_OpAmp_Vminus") # FET Source va all'ingresso invertente


    def _connect_nodes(self):
        """Connette i nodi dell'intero circuito."""
        # Connessione dell'Audio Input al VCA
        # (Già fatto in _add_components per R_in)

        # L'uscita dell'Op-Amp VCA è l'uscita audio del compressore
        self.connect_nodes("Audio_Output", "VCA_OpAmp_Vout")

        # Connessione dell'audio input del compressore all'input della sidechain
        self.connect_nodes("Audio_Input", "SC_Audio_Input")

        # La tensione di controllo dalla side-chain va al Gate del FET
        self.connect_nodes("SC_Control_Voltage_Output", "FET_Gate")
        
        # Il source del FET è già connesso all'OpAmp_Vminus, e l'OpAmp_Vminus è a virtual GND.
        # Questo è importante per il bias del FET.

    def _setup_sidechain_times(self, attack_ms, release_ms):
        """
        Imposta i tempi di attack e release della side-chain.
        Questo metodo deve essere chiamato dopo che il circuito è stato completamente costruito
        e i valori R e C della side-chain sono stati impostati nel sottocircuito.
        """
        # Calcoliamo i valori di R_attack e R_load per la side-chain in base ai tempi desiderati.
        # Assumiamo un valore fisso per il condensatore di release della side-chain per semplicità.
        # Questo valore è anche quello usato nel DiodeBridgeCircuit interno.
        C_release_fixed = 1.0e-6 # Un microfarad come valore comune

        R_attack_calc = (attack_ms / 1000.0) / C_release_fixed
        R_load_calc = (release_ms / 1000.0) / C_release_fixed

        # Aggiorniamo i parametri della sidechain e l'istanza.
        self.sidechain_params["attack_resistor"] = R_attack_calc
        self.sidechain_params["release_capacitor"] = C_release_fixed
        self.sidechain_params["load_resistor"] = R_load_calc

        # Ora creiamo l'istanza finale della side-chain con i parametri calcolati
        self.sidechain = DiodeBridgeCompressorSidechain(
            diode_is=self.sidechain_params["diode_is"],
            diode_n=self.sidechain_params["diode_n"],
            attack_resistor=self.sidechain_params["attack_resistor"],
            release_capacitor=self.sidechain_params["release_capacitor"],
            load_resistor=self.sidechain_params["load_resistor"],
            sample_rate=self.sample_rate
        )

        print(f"Side-chain: Attack R={self.sidechain_params['attack_resistor']:.1f} Ohm, "
              f"Release C={self.sidechain_params['release_capacitor']:.1e} Farad, "
              f"Load R={self.sidechain_params['load_resistor']:.1f} Ohm")


    def get_audio_input_node(self):
        return self.get_node_id("Audio_Input")

    def get_audio_output_node(self):
        return self.get_node_id("Audio_Output")

    def get_sidechain_input_node(self):
        return self.get_node_id("SC_Audio_Input")

    def get_sidechain_control_voltage_output_node(self):
        return self.get_node_id("SC_Control_Voltage_Output")

    def get_fet_instance(self):
        return self.fet

    def get_sidechain_instance(self):
        return self.sidechain

    def set_attack_release_times(self, attack_ms, release_ms):
        """Metodo esterno per impostare i tempi di attack e release della side-chain."""
        self._setup_sidechain_times(attack_ms, release_ms)
        # Qui potresti voler ri-inizializzare il circuito o informare il solver.
        # Per ora, i componenti interni del solver (nel main loop) devono essere aggiornati.
