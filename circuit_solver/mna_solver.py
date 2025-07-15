# circuit_solver/mna_solver.py

import numpy as np
from scipy.optimize import fsolve
from components.component import Component
from components.voltage_source import VoltageSource
from components.diode import Diode
from components.mosfet import MOSFET
from components.triode import Triode
from components.pentode import Pentode
from components.rectifier_tube import RectifierTube
from components.led import LED
from components.ldr import LDR
from components.jfet import JFET
from components.bjt import BJT
from components.speaker_driver import SpeakerDriver
from components.closed_box_cabinet import ClosedBoxCabinet
from components.bass_reflex_cabinet import BassReflexCabinet
from components.delay_line import DelayLine
from components.splitter import Splitter # Importa il nuovo Splitter

class MnaSolver:
    """
    Risolve il circuito usando l'Analisi Nodal Modificata (MNA) per simulazioni transitorie.
    """
    def __init__(self, circuit):
        self.circuit = circuit
        self.num_nodes = circuit.get_num_nodes()
        self.num_voltage_sources = len(circuit.get_voltage_sources())
        self.num_splitter_aux_vars = sum(s.num_outputs for s in circuit.get_splitters()) # Variabili ausiliarie per splitter
        
        # Dimensione totale del sistema di equazioni MNA
        self.num_total_equations = self.num_nodes + self.num_voltage_sources + self.num_splitter_aux_vars

        # Assegna gli indici delle variabili ausiliarie alle sorgenti di tensione e agli splitter
        self._assign_auxiliary_indices()

        # Classifica i componenti per un accesso più rapido nel solutore
        self.linear_components = []
        self.dynamic_components = [] # Componenti con stato (C, L)
        self.nonlinear_components = [] # Componenti non lineari (Diodo, MOSFET, Triodo, ecc.)
        self.functional_components = [] # Componenti funzionali (DelayLine, LDR, ecc.)

        self._classify_components()

    def _assign_auxiliary_indices(self):
        """
        Assegna gli indici nel vettore delle incognite alle correnti delle sorgenti di tensione
        e alle correnti ausiliarie degli splitter.
        """
        current_aux_idx = self.num_nodes # Inizia dopo gli ID dei nodi

        # Assegna indici alle correnti delle sorgenti di tensione
        for vs in self.circuit.get_voltage_sources():
            vs._set_current_index(current_aux_idx)
            current_aux_idx += 1
        
        # Assegna indici alle correnti ausiliarie degli splitter
        for splitter in self.circuit.get_splitters():
            indices = []
            for _ in range(splitter.num_outputs):
                indices.append(current_aux_idx)
                current_aux_idx += 1
            splitter._set_output_current_indices(indices)

    def _classify_components(self):
        """Classifica i componenti in base al loro tipo per una gestione efficiente."""
        for comp in self.circuit.get_components():
            if isinstance(comp, (Diode, MOSFET, Triode, Pentode, RectifierTube, LED, JFET, BJT)):
                self.nonlinear_components.append(comp)
            elif isinstance(comp, (DelayLine, LDR)): # LDR è dinamico ma il suo stamp è lineare
                self.functional_components.append(comp)
            elif hasattr(comp, 'update_state'): # Condensatori, Induttori, SpeakerDriver, Cabinets
                self.dynamic_components.append(comp)
            else: # Resistori, Sorgenti di Corrente, Sorgenti di Tensione, Splitter
                self.linear_components.append(comp)

    def _system_equations(self, x: np.ndarray, dt: float, prev_solution: np.ndarray, time: float) -> np.ndarray:
        """
        Definisce il sistema di equazioni MNA (F(x) = 0) per fsolve.
        Args:
            x (np.ndarray): Il vettore delle incognite (tensioni ai nodi, correnti delle Vs, correnti ausiliarie).
            dt (float): Il passo temporale.
            prev_solution (np.ndarray): La soluzione del passo temporale precedente.
            time (float): Il tempo attuale della simulazione.
        Returns:
            np.ndarray: Il vettore delle equazioni.
        """
        # Inizializza la matrice MNA (A) e il vettore RHS (B)
        # Questi saranno usati per costruire la parte lineare del sistema
        A = np.zeros((self.num_total_equations, self.num_total_equations))
        B = np.zeros(self.num_total_equations)

        # --- Contributi dei componenti lineari e dinamici ---
        # Questi vengono aggiunti alla matrice A e al vettore B
        for comp in self.linear_components + self.dynamic_components:
            stamp_A, stamp_B = comp.get_stamps(self.num_total_equations, dt, x, prev_solution, time)
            A += stamp_A
            B += stamp_B

        # --- Contributi dei componenti non lineari (gestiti implicitamente) ---
        # Per i componenti non lineari, le loro equazioni sono aggiunte direttamente qui.
        # F(x) = A*x - B - I_nonlinear(x) = 0
        # Quindi, I_nonlinear(x) viene sottratta dal vettore B.
        
        # Inizializza il vettore delle correnti non lineari
        I_nonlinear = np.zeros(self.num_total_equations)

        for comp in self.nonlinear_components:
            node_ids = comp.node_ids
            # Recupera le tensioni ai capi del componente non lineare
            if isinstance(comp, Diode) or isinstance(comp, SchottkyDiode) or isinstance(comp, ZenerDiode):
                Vd = x[node_ids[0]] - x[node_ids[1]]
                current = comp.calculate_current(Vd)
                if node_ids[0] != 0: I_nonlinear[node_ids[0]] += current
                if node_ids[1] != 0: I_nonlinear[node_ids[1]] -= current
            elif isinstance(comp, MOSFET):
                Vds = x[node_ids[0]] - x[node_ids[2]] # Drain - Source
                Vgs = x[node_ids[1]] - x[node_ids[2]] # Gate - Source
                Id = comp.calculate_drain_current(Vgs, Vds)
                if node_ids[0] != 0: I_nonlinear[node_ids[0]] += Id # Corrente entra nel Drain
                if node_ids[2] != 0: I_nonlinear[node_ids[2]] -= Id # Corrente esce dal Source
            elif isinstance(comp, Triode):
                Vpk = x[node_ids[0]] - x[node_ids[2]] # Plate - Cathode
                Vgk = x[node_ids[1]] - x[node_ids[2]] # Grid - Cathode
                Ip = comp.calculate_plate_current(Vgk, Vpk)
                if node_ids[0] != 0: I_nonlinear[node_ids[0]] += Ip # Corrente entra nella Placca
                if node_ids[2] != 0: I_nonlinear[node_ids[2]] -= Ip # Corrente esce dal Catodo
            elif isinstance(comp, Pentode):
                Vpk = x[comp.node_map['plate']] - x[comp.node_map['cathode']]
                Vgk = x[comp.node_map['grid']] - x[comp.node_map['cathode']]
                Vg2k = x[comp.node_map['screen_grid']] - x[comp.node_map['cathode']]
                Vg3k = x[comp.node_map['suppressor_grid']] - x[comp.node_map['cathode']]
                Ip = comp.calculate_plate_current(Vgk, Vpk, Vg2k, Vg3k)
                if comp.node_map['plate'] != 0: I_nonlinear[comp.node_map['plate']] += Ip
                if comp.node_map['cathode'] != 0: I_nonlinear[comp.node_map['cathode']] -= Ip
            elif isinstance(comp, RectifierTube):
                Vd = x[node_ids[0]] - x[node_ids[1]] # Anode - Cathode
                current = comp.calculate_current(Vd)
                if node_ids[0] != 0: I_nonlinear[node_ids[0]] += current
                if node_ids[1] != 0: I_nonlinear[node_ids[1]] -= current
            elif isinstance(comp, BJT):
                Vbe = x[node_ids[1]] - x[node_ids[2]] # Base - Emitter
                Vbc = x[node_ids[1]] - x[node_ids[0]] # Base - Collector (nota: Vbc, non Vcb)
                Ic, Ib = comp.calculate_currents(Vbe, Vbc)
                Ie = Ic + Ib # Corrente di emettitore
                
                if comp.type == 'npn':
                    if node_ids[0] != 0: I_nonlinear[node_ids[0]] += Ic # Ic entra nel Collector
                    if node_ids[1] != 0: I_nonlinear[node_ids[1]] += Ib # Ib entra nella Base
                    if node_ids[2] != 0: I_nonlinear[node_ids[2]] -= Ie # Ie esce dall'Emitter
                elif comp.type == 'pnp':
                    if node_ids[0] != 0: I_nonlinear[node_ids[0]] -= Ic # Ic esce dal Collector
                    if node_ids[1] != 0: I_nonlinear[node_ids[1]] -= Ib # Ib esce dalla Base
                    if node_ids[2] != 0: I_nonlinear[node_ids[2]] += Ie # Ie entra nell'Emitter
            elif isinstance(comp, JFET):
                Vds = x[node_ids[0]] - x[node_ids[2]] # Drain - Source
                Vgs = x[node_ids[1]] - x[node_ids[2]] # Gate - Source
                Id = comp.calculate_drain_current(Vgs, Vds)
                if node_ids[0] != 0: I_nonlinear[node_ids[0]] += Id # Id entra nel Drain
                if node_ids[2] != 0: I_nonlinear[node_ids[2]] -= Id # Id esce dal Source
            # LDR è gestito come un resistore dinamico, quindi è in linear_components.
            # OpAmp è gestito come sorgente di tensione controllata.
            # SpeakerDriver, Cabinets sono in dynamic_components.

        # Costruisci l'equazione finale: A*x - B_total = 0
        # B_total = B_linear + I_nonlinear (correnti che escono dai nodi, quindi sottratte)
        # Equazioni MNA: G*V + C*dV/dt + I_nl = I_s
        # Con integrazione trapezoidale: G*V_n + (2C/dt)*V_n - (2C/dt)*V_n-1 - I_n-1 = I_s
        # F(x) = A*x - B_linear - I_nonlinear = 0
        
        # Il vettore B contiene già i contributi delle sorgenti indipendenti e dei termini dinamici.
        # Ora sottraiamo le correnti non lineari (che sono funzioni di x)
        
        # Per fsolve, vogliamo F(x) = 0.
        # Il sistema lineare è A*x = B.
        # Per i componenti non lineari, la loro corrente è una funzione di x.
        # Quindi, l'equazione KCL per un nodo diventa: sum(correnti_lineari) + I_non_lineare(x) = 0
        # Nel nostro sistema A*x = B, questo si traduce in:
        # A_linear * x - B_linear - I_non_lineare(x) = 0
        # Quindi, il vettore di residui è (A @ x) - B - I_nonlinear
        
        # Nota: il nodo 0 (ground) non ha un'equazione KCL.
        # La riga 0 della matrice A e del vettore B non viene utilizzata per KCL.
        # Per convenzione, V_ground = 0.
        
        # La riga 0 di F è sempre 0, per mantenere la dimensione.
        F = (A @ x) - B - I_nonlinear
        F[0] = x[0] # Impone V_ground = 0

        return F

    def simulate_transient(self, start_time: float, end_time: float, time_step: float) -> tuple:
        """
        Esegue una simulazione transitoria del circuito.
        Args:
            start_time (float): Tempo di inizio della simulazione.
            end_time (float): Tempo di fine della simulazione.
            time_step (float): Passo temporale della simulazione (dt).
        Returns:
            tuple: (times, solution_history) - Array dei tempi e matrice delle soluzioni.
        """
        times = np.arange(start_time, end_time + time_step, time_step)
        
        # Inizializza la soluzione precedente a zero per tutte le incognite
        prev_solution = np.zeros(self.num_total_equations)
        solution_history = []

        print(f"Inizio simulazione transitoria da {start_time}s a {end_time}s con dt={time_step}s.")
        print(f"Numero totale di equazioni MNA: {self.num_total_equations}")

        for i, t in enumerate(times):
            # Guess iniziale per fsolve (soluzione del passo precedente)
            initial_guess = prev_solution.copy()

            # Risolvi il sistema di equazioni non lineari per il passo attuale
            try:
                # Passa prev_solution e time a _system_equations
                current_solution = fsolve(self._system_equations, initial_guess, args=(time_step, prev_solution, t))
                
                # Aggiorna lo stato dei componenti dinamici per il prossimo passo
                self._update_dynamic_component_states(current_solution, prev_solution, time_step)
                
                prev_solution = current_solution
                solution_history.append(current_solution)
                
                # print(f"Tempo: {t:.6f}s, Soluzione: {current_solution[self.circuit.node_map.get('output_final', 0)]:.4f}V")

            except Exception as e:
                print(f"Errore durante la risoluzione a t={t:.6f}s: {e}")
                solution_history.append(np.full(self.num_total_equations, np.nan)) # Aggiungi NaN in caso di errore
                break
        
        return np.array(times[:len(solution_history)]), np.array(solution_history)

    def _update_dynamic_component_states(self, current_solution: np.ndarray, prev_solution: np.ndarray, dt: float):
        """
        Aggiorna lo stato interno dei componenti dinamici (C, L, SpeakerDriver, Cabinets).
        """
        for comp in self.dynamic_components:
            node1_id, node2_id = comp.node_ids[0], comp.node_ids[1] # Assumi i primi due nodi per V e I

            V_curr = current_solution[node1_id] - current_solution[node2_id]
            V_prev = prev_solution[node1_id] - prev_solution[node2_id]

            if isinstance(comp, (Capacitor, ClosedBoxCabinet, BassReflexCabinet)):
                # Corrente attraverso il condensatore (trapezoidale)
                I_curr = (2.0 * comp.capacitance / dt) * (V_curr - V_prev) - comp.i_prev
                comp.update_state(V_curr, I_curr)
            elif isinstance(comp, Inductor):
                # Corrente attraverso l'induttore (trapezoidale)
                I_curr = comp.i_prev + (dt / (2.0 * comp.L)) * (V_curr + V_prev)
                comp.update_state(V_curr, I_curr)
            elif isinstance(comp, SpeakerDriver):
                # Per SpeakerDriver, i nodi sono (elec_plus, elec_minus, mech_velocity_out)
                # Le, L_mech, C_mech hanno i loro stati.
                # Questa parte è più complessa e richiederebbe nodi interni espliciti per essere accurata.
                # Per ora, aggiorniamo solo i parametri che sono stati esplicitamente modellati come dinamici nel get_stamps.
                # Se i get_stamps di SpeakerDriver sono stati semplificati, l'update_state potrebbe essere vuoto o parziale.
                pass # L'update_state di SpeakerDriver è più complesso e dipende dai suoi nodi interni.

        # Aggiorna lo stato dei componenti funzionali (es. DelayLine)
        # Questi non sono gestiti da MNA ma dalla pipeline.
        # Se un DelayLine fosse direttamente nel circuito MNA (come sorgente V/I controllata)
        # allora il suo update sarebbe qui.
        # Per ora, DelayLine è gestito esternamente nel main.
        # LDR è gestito come resistore dinamico, il suo stato è la resistenza attuale.
        pass

