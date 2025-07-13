# circuit_solver/mna_solver.py

import numpy as np
from scipy.optimize import fsolve

# Assicurati che queste classi siano importate correttamente o definite nello stesso file
# Se i tuoi componenti sono in 'components/', assicurati che il percorso sia corretto.
from components.resistor import Resistor
from components.capacitor import Capacitor
# from components.inductor import Inductor # Se hai implementato induttori
from components.voltage_source import VoltageSource
from components.current_source import CurrentSource
from components.diode import Diode
from components.mosfet import MOSFET
from components.triode import Triode
from components.op_amp import OpAmp # La tua classe OpAmp (nota il nome del file op_amp.py)

# --- Nuova Classe: DelayLine (Componente Funzionale) ---
# Questa classe è inclusa qui per comodità, ma potrebbe essere in un file separato (es. components/delay_line.py)
class DelayLine:
    """
    Simula una linea di ritardo discreta in termini di campioni.
    Questa è una componente funzionale, non un componente SPICE tradizionale.
    """
    def __init__(self, max_delay_seconds, sample_rate):
        self.sample_rate = sample_rate
        self.max_samples = int(max_delay_seconds * sample_rate)
        if self.max_samples < 1:
            self.max_samples = 1 # Minimum 1 sample delay
        self.buffer = np.zeros(self.max_samples) # Inizializza con zeri
        self.write_index = 0
        self.current_delay_samples = 0 # Ritardo attuale in campioni
        
        # Storico dell'input e output per debugging o uso interno
        self._input_history = []
        self._output_history = []

    def set_delay_time(self, delay_time_seconds):
        """Imposta il tempo di ritardo in secondi."""
        num_samples = int(delay_time_seconds * self.sample_rate)
        if num_samples > self.max_samples:
            num_samples = self.max_samples
        self.current_delay_samples = num_samples
    
    def get_delay_time(self):
        return self.current_delay_samples / self.sample_rate

    def update(self, input_sample):
        """
        Aggiorna la linea di ritardo con un nuovo campione e restituisce il campione ritardato.
        """
        # Scrivi il nuovo campione nel buffer
        self.buffer[self.write_index] = input_sample
        self._input_history.append(input_sample) # Per debugging

        # Calcola l'indice di lettura
        read_index = (self.write_index - self.current_delay_samples + self.max_samples) % self.max_samples
        
        # Ottieni il campione ritardato
        output_sample = self.buffer[read_index]
        self._output_history.append(output_sample) # Per debugging

        # Avanza l'indice di scrittura (buffer circolare)
        self.write_index = (self.write_index + 1) % self.max_samples
        
        return output_sample

    def reset(self):
        """Resetta il buffer del delay a zero."""
        self.buffer = np.zeros(self.max_samples)
        self.write_index = 0
        self._input_history = []
        self._output_history = []

# --- MnaSolver aggiornato ---
class MnaSolver:
    def __init__(self, circuit):
        self.circuit = circuit
        self._node_map = circuit.get_node_map()
        self._num_nodes = circuit.get_num_nodes()
        self._num_voltage_sources = circuit.get_num_voltage_sources()
        
        # Gli OpAmp e altri componenti non lineari sono gestiti da fsolve.
        # Non aggiungono equazioni ausiliarie lineari extra al conteggio principale di _total_equations.
        self._num_nonlinear_components = len(circuit.nonlinear_components)
        
        # Il numero totale di equazioni/incognite include solo nodi e correnti delle sorgenti di tensione.
        self._total_equations = self._num_nodes + self._num_voltage_sources

        # Mappatura per le sorgenti di tensione (necessaria per le loro colonne/righe in MNA)
        self._voltage_source_map = {}
        vs_idx = 0
        for i, comp in enumerate(circuit.components):
            if isinstance(comp, VoltageSource):
                self._voltage_source_map[comp.name] = vs_idx
                vs_idx += 1
        
        # Riferimento alla DelayLine (se presente nel circuito)
        self.delay_line_instance = None
        if hasattr(circuit, 'delay_line'):
            self.delay_line_instance = circuit.delay_line
        
        self.initial_conditions = None

    def _get_node_index(self, node_name):
        return self._node_map.get(node_name, -1) # Restituisce -1 se il nodo non è trovato (es. 'gnd')

    def _get_voltage_source_index(self, vs_name):
        return self._voltage_source_map.get(vs_name, -1)

    def _build_linear_system(self, x_current_guess, is_transient=False, dt=0, prev_x=None, time=0):
        """
        Costruisce il sistema lineare MNA A*x = b per i componenti lineari e dinamici.
        I componenti non lineari (Diode, MOSFET, Triode, OpAmp) sono gestiti in _system_equations.
        x_current_guess: il guess corrente per le tensioni ai nodi e le correnti delle sorgenti di tensione.
                         Usato per calcolare i termini dinamici dei condensatori.
        """
        
        A_matrix = np.zeros((self._total_equations, self._total_equations))
        b_vector = np.zeros(self._total_equations)

        # Offset per le sezioni della matrice/vettore relative alle sorgenti di tensione
        vs_offset = self._num_nodes 
        
        for comp in self.circuit.components:
            n1_idx = self._get_node_index(comp.node1) if hasattr(comp, 'node1') else -1
            n2_idx = self._get_node_index(comp.node2) if hasattr(comp, 'node2') else -1

            if isinstance(comp, Resistor):
                g = 1.0 / comp.resistance
                if n1_idx != -1: A_matrix[n1_idx, n1_idx] += g
                if n2_idx != -1: A_matrix[n2_idx, n2_idx] += g
                if n1_idx != -1 and n2_idx != -1:
                    A_matrix[n1_idx, n2_idx] -= g
                    A_matrix[n2_idx, n1_idx] -= g

            elif isinstance(comp, VoltageSource):
                vs_col_idx = self._get_voltage_source_index(comp.name) + vs_offset # Colonna per la corrente Vs
                vs_row_idx = self._get_voltage_source_index(comp.name) + vs_offset # Riga per l'equazione Vs
                
                # Contributo alla KCL ai nodi (corrente della Vs)
                if n1_idx != -1: A_matrix[n1_idx, vs_col_idx] += 1
                if n2_idx != -1: A_matrix[n2_idx, vs_col_idx] -= 1

                # Equazione della sorgente di tensione: V_n1 - V_n2 = V_source
                if n1_idx != -1: A_matrix[vs_row_idx, n1_idx] += 1
                if n2_idx != -1: A_matrix[vs_row_idx, n2_idx] -= 1
                b_vector[vs_row_idx] = comp.get_voltage(time)

            elif isinstance(comp, CurrentSource):
                # Aggiungi la corrente al vettore b
                if n1_idx != -1: b_vector[n1_idx] -= comp.current # Corrente che esce dal nodo n1
                if n2_idx != -1: b_vector[n2_idx] += comp.current # Corrente che entra nel nodo n2

            elif is_transient and isinstance(comp, Capacitor):
                # Modello per l'integrazione trapezoidale implicita di un condensatore
                g_eq = 2.0 * comp.capacitance / dt
                
                # Aggiungi la conduttanza equivalente
                if n1_idx != -1: A_matrix[n1_idx, n1_idx] += g_eq
                if n2_idx != -1: A_matrix[n2_idx, n2_idx] += g_eq
                if n1_idx != -1 and n2_idx != -1:
                    A_matrix[n1_idx, n2_idx] -= g_eq
                    A_matrix[n2_idx, n1_idx] -= g_eq

                # Calcola V_C(t-dt) = V_n1(t-dt) - V_n2(t-dt)
                V_C_prev = 0.0
                if n1_idx != -1 and prev_x is not None: V_C_prev += prev_x[n1_idx]
                if n2_idx != -1 and prev_x is not None: V_C_prev -= prev_x[n2_idx]

                # Aggiungi la sorgente di corrente equivalente al vettore b
                i_eq = g_eq * V_C_prev
                
                if n1_idx != -1: b_vector[n1_idx] -= i_eq # Corrente che esce da n1
                if n2_idx != -1: b_vector[n2_idx] += i_eq # Corrente che entra in n2
            
            # I componenti non lineari (OpAmp, Diode, MOSFET, Triode) non vengono aggiunti qui,
            # ma la loro logica è in _system_equations.
            # La resistenza di ingresso dell'OpAmp (Rin) dovrebbe essere modellata come un Resistor separato
            # nel circuito che usa l'OpAmp, tra i nodi in_inv e in_non_inv.

        return A_matrix, b_vector

    def _system_equations(self, x, A_linear, b_linear, dt, prev_x, time):
        """
        Sistema di equazioni MNA, inclusi i contributi dei componenti non lineari.
        x: vettore delle incognite (tensioni ai nodi + correnti sorgenti di tensione)
        A_linear, b_linear: matrici MNA per i componenti lineari
        """
        
        # Calcola il residuo iniziale del sistema lineare: A*x - b
        F = np.dot(A_linear, x) - b_linear
        
        # Aggiungi i contributi delle non linearità al vettore F
        # Ogni componente non lineare contribuisce con correnti ai nodi.
        # Per fsolve, vogliamo F(x) = 0.
        
        for comp in self.circuit.nonlinear_components:
            # Assicurati che i nodi siano mappati correttamente nel circuito padre.
            n1_idx = self._get_node_index(comp.node1) if hasattr(comp, 'node1') else -1
            n2_idx = self._get_node_index(comp.node2) if hasattr(comp, 'node2') else -1

            if isinstance(comp, Diode):
                V_anode = x[n1_idx] if n1_idx != -1 else 0.0
                V_cathode = x[n2_idx] if n2_idx != -1 else 0.0
                
                I_diode = comp.calculate_current(V_anode - V_cathode)
                
                if n1_idx != -1: F[n1_idx] -= I_diode # Corrente che esce dall'anodo
                if n2_idx != -1: F[n2_idx] += I_diode # Corrente che entra nel catodo

            elif isinstance(comp, MOSFET):
                V_drain = x[self._get_node_index(comp.nodes['drain'])] if 'drain' in comp.nodes else 0.0
                V_gate = x[self._get_node_index(comp.nodes['gate'])] if 'gate' in comp.nodes else 0.0
                V_source = x[self._get_node_index(comp.nodes['source'])] if 'source' in comp.nodes else 0.0
                
                Vgs = V_gate - V_source
                Vds = V_drain - V_source
                
                Id = comp.calculate_drain_current(Vgs, Vds)
                
                drain_idx = self._get_node_index(comp.nodes['drain'])
                source_idx = self._get_node_index(comp.nodes['source'])

                if drain_idx != -1: F[drain_idx] -= Id # Corrente che esce dal drain
                if source_idx != -1: F[source_idx] += Id # Corrente che entra nel source

            elif isinstance(comp, Triode):
                V_plate = x[self._get_node_index(comp.nodes['plate'])] if 'plate' in comp.nodes else 0.0
                V_grid = x[self._get_node_index(comp.nodes['grid'])] if 'grid' in comp.nodes else 0.0
                V_cathode = x[self._get_node_index(comp.nodes['cathode'])] if 'cathode' in comp.nodes else 0.0

                Vgk = V_grid - V_cathode
                Vpk = V_plate - V_cathode

                Ip = comp.calculate_plate_current(Vgk, Vpk)

                plate_idx = self._get_node_index(comp.nodes['plate'])
                grid_idx = self._get_node_index(comp.nodes['grid'])
                cathode_idx = self._get_node_index(comp.nodes['cathode'])

                if plate_idx != -1: F[plate_idx] -= Ip # Corrente che esce dalla placca
                if cathode_idx != -1: F[cathode_idx] += Ip # Corrente che entra nel catodo
            
            elif isinstance(comp, OpAmp):
                # La tua classe OpAmp è un componente non lineare a causa della saturazione.
                # Deve essere trattato come una sorgente di corrente che dipende dalle tensioni dei nodi.
                
                in_inv_idx = self._get_node_index(comp.nodes['in_inv'])
                in_non_inv_idx = self._get_node_index(comp.nodes['in_non_inv'])
                out_idx = self._get_node_index(comp.nodes['out'])
                
                V_in_inv = x[in_inv_idx] if in_inv_idx != -1 else 0.0
                V_in_non_inv = x[in_non_inv_idx] if in_non_inv_idx != -1 else 0.0
                V_out_circuit = x[out_idx] if out_idx != -1 else 0.0 # Tensione del nodo di uscita

                # Calcola la tensione di uscita ideale (con saturazione) dell'OpAmp
                V_out_opamp_internal = comp.calculate_output_voltage(V_in_non_inv, V_in_inv)

                # La corrente che l'OpAmp inietta nel nodo di uscita 'out_idx' è data da:
                # I_out = (V_out_opamp_internal - V_out_circuit) / Rout
                # Questa corrente ESCE dal nodo di output dell'OpAmp.
                # Quindi, contribuisce negativamente alla KCL del nodo di output.
                if out_idx != -1:
                    # F[out_idx] rappresenta la somma delle correnti che ESCONO dal nodo.
                    # Quindi, la corrente che l'OpAmp inietta (che esce dal suo pin di output)
                    # deve essere SOTTRATTA da F[out_idx] per mantenere F[out_idx] = 0.
                    # Se consideriamo la corrente che ENTRA nel nodo, allora sarebbe aggiunta.
                    # Per convenzione MNA, I_component = (V_node1 - V_node2) / R.
                    # Se la corrente esce da node1, è positiva.
                    # Qui, la corrente esce dall'OpAmp verso il nodo out_idx.
                    # Quindi, nel residuo F[out_idx] (che è KCL = Sum of currents leaving the node),
                    # dobbiamo sottrarre questa corrente.
                    F[out_idx] -= (V_out_opamp_internal - V_out_circuit) / comp.Rout
                
                # Le correnti di input (attraverso Rin) sono gestite come un resistore separato nel circuito.
                # Quindi non c'è bisogno di aggiungere termini a F[in_inv_idx] o F[in_non_inv_idx] qui.

        return F

    def solve_dc(self, initial_guess=None):
        """
        Risolve il circuito in DC (stato stazionario).
        Per i componenti dinamici (C, L): C = circuito aperto, L = corto circuito.
        """
        
        # Inizializza x_guess
        if initial_guess is None:
            initial_guess = np.zeros(self._total_equations)
        else:
            if len(initial_guess) != self._total_equations:
                 print(f"Warning: Initial guess size mismatch. Expected {self._total_equations}, got {len(initial_guess)}. Resetting to zeros.")
                 initial_guess = np.zeros(self._total_equations)

        # Build linear part of MNA (independent of non-linearities) for DC analysis
        A_linear, b_linear = self._build_linear_system(initial_guess, is_transient=False, dt=0, prev_x=None, time=0)
        
        try:
            # Wrap _system_equations to pass fixed arguments
            def func(x):
                # Per DC, dt=0 e prev_x=None
                return self._system_equations(x, A_linear, b_linear, dt=0, prev_x=None, time=0)

            # Solve the nonlinear system using fsolve
            solution, info, ier, msg = fsolve(func, initial_guess, full_output=True)

            if ier == 1:
                # Extract node voltages and voltage source currents
                node_voltages = solution[:self._num_nodes]
                source_currents = solution[self._num_nodes:self._num_nodes + self._num_voltage_sources]
                return solution, node_voltages, source_currents # Restituisce la soluzione completa per transient
            else:
                print(f"DC Convergence failed: {msg}")
                return None, None, None
        except Exception as e:
            print(f"Error during DC simulation: {e}")
            return None, None, None

    def solve_transient(self, t_start, t_end, dt, output_nodes=None):
        """
        Esegue la simulazione transitoria del circuito.
        t_start: tempo di inizio simulazione
        t_end: tempo di fine simulazione
        dt: passo temporale
        output_nodes: lista di nomi di nodi da monitorare (se None, monitora tutti i nodi)
        """
        print(f"Starting transient simulation from {t_start}s to {t_end}s with dt={dt}s...")

        time_points = np.arange(t_start, t_end + dt, dt)
        
        # Inizializzazione della soluzione al tempo t-dt (condizioni iniziali)
        # Prima, calcola il punto di riposo DC come condizione iniziale.
        prev_x_solution_full, _, _ = self.solve_dc()
        if prev_x_solution_full is None:
            print("Failed to converge DC solution for initial conditions. Aborting transient.")
            return None, None

        x_prev = np.copy(prev_x_solution_full) # Usa la soluzione DC completa come guess iniziale per il primo passo

        # Vettore per salvare i risultati
        results = {node_name: [] for node_name in (output_nodes if output_nodes else self._node_map.keys())}
        results['time'] = []
        
        # Assicurati che la DelayLine sia resettata prima della simulazione
        if self.delay_line_instance:
            self.delay_line_instance.reset()
            # Imposta il tempo di ritardo iniziale, se applicabile al circuito
            # self.delay_line_instance.set_delay_time(self.circuit.delay_time_seconds) # Esempio

        for i, t in enumerate(time_points):
            # Aggiorna i valori delle sorgenti di tensione e corrente dipendenti dal tempo
            for comp in self.circuit.components:
                if isinstance(comp, VoltageSource):
                    comp.set_voltage(comp.get_voltage(t)) # La get_voltage() deve essere implementata per dipendere dal tempo
                if isinstance(comp, CurrentSource):
                    comp.set_current(comp.get_current(t)) # Simile per la corrente

            # Se hai una DelayLine, devi aggiornare la sua uscita nel circuito
            # Questo è un punto di integrazione cruciale.
            # Assumiamo che il circuito abbia un nodo 'delay_input' e una sorgente 'V_delay_output_source'
            if self.delay_line_instance:
                # Preleva la tensione del nodo di input del delay dal passo precedente
                delay_input_node_idx = self._get_node_index('delay_input')
                if delay_input_node_idx != -1:
                    delay_input_voltage = x_prev[delay_input_node_idx] 
                    delayed_voltage = self.delay_line_instance.update(delay_input_voltage)
                    
                    # Ora, aggiorna la sorgente di tensione "virtuale" che rappresenta l'uscita del delay
                    # Il nome della sorgente deve corrispondere a quella nel circuito
                    delay_output_source_comp = next((c for c in self.circuit.components if c.name == "V_delay_output_source"), None)
                    if delay_output_source_comp:
                        delay_output_source_comp.set_voltage(delayed_voltage) # Imposta la tensione della sorgente
                    else:
                        print(f"Warning: DelayLine output source 'V_delay_output_source' not found in circuit at t={t:.6f}s.")
                else:
                    print(f"Warning: DelayLine input node 'delay_input' not found in circuit at t={t:.6f}s.")


            # Build linear part of MNA for this timestep
            A_linear, b_linear = self._build_linear_system(x_prev, is_transient=True, dt=dt, prev_x=x_prev, time=t)

            # Define the function for fsolve for the current timestep
            def func_timestep(x_current):
                return self._system_equations(x_current, A_linear, b_linear, dt=dt, prev_x=x_prev, time=t)

            try:
                # Usa x_prev come guess iniziale per il passo corrente
                x_current_solution, info, ier, msg = fsolve(func_timestep, x_prev, full_output=True)

                if ier == 1:
                    x_prev = x_current_solution # Aggiorna la soluzione per il prossimo passo
                    
                    # Salva i risultati per i nodi specificati
                    results['time'].append(t)
                    for node_name in (output_nodes if output_nodes else self._node_map.keys()):
                        node_idx = self._get_node_index(node_name)
                        if node_idx != -1:
                            results[node_name].append(x_current_solution[node_idx])
                        else:
                            results[node_name].append(np.nan) # Nodo non trovato

                else:
                    print(f"Transient simulation failed to converge at t={t:.6f}s: {msg}")
                    break # Interrompi la simulazione in caso di non convergenza

            except Exception as e:
                print(f"Error during transient simulation at t={t:.6f}s: {e}")
                break

        print("Transient simulation finished.")
        return results, time_points
