# circuit_solver/mna_solver.py

import numpy as np
from scipy.optimize import fsolve # Per la risoluzione di sistemi non lineari

# Assicurati che queste classi siano importate correttamente o definite nello stesso file
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor # Se hai implementato induttori
from components.voltage_source import VoltageSource
from components.current_source import CurrentSource # Se hai implementato sorgenti di corrente
from components.diode import Diode
from components.mosfet import MOSFET
from components.triode import Triode
from components.opamp import OpAmp # Nuova classe OpAmp


# --- Nuova Classe: DelayLine (Componente Funzionale) ---
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
        self.current_delay_samples = 0
        
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
        self._num_nonlinear_components = len(circuit.nonlinear_components)
        self._total_equations = self._num_nodes + self._num_voltage_sources

        # Mappatura per le sorgenti di tensione per MNA (per tracciare quali equazioni sono associate a Vs)
        self._voltage_source_map = {}
        vs_idx = 0
        for i, comp in enumerate(circuit.components):
            if isinstance(comp, VoltageSource):
                self._voltage_source_map[comp.name] = vs_idx
                vs_idx += 1
        
        # Riferimento alla DelayLine se presente nel circuito
        self.delay_line_instance = None
        # Trova la DelayLine nel circuito se è stata aggiunta come attributo
        # Questa è una convenzione, potresti voler gestire le DelayLine in modo più generico
        if hasattr(circuit, 'delay_line'): # Assumi che il circuito abbia un attributo 'delay_line'
            self.delay_line_instance = circuit.delay_line
        
        # Per la simulazione transitoria
        self.initial_conditions = None

    def _get_node_index(self, node_name):
        return self._node_map.get(node_name, -1) # Restituisce -1 se il nodo non è trovato (es. 'gnd')

    def _get_voltage_source_index(self, vs_name):
        return self._voltage_source_map.get(vs_name, -1)

    def _build_linear_system(self, x, G, B, C, D, is_transient=False, dt=0, prev_x=None, time=0):
        """
        Costruisce il sistema lineare MNA A*x = b.
        x: vettore delle incognite (tensioni ai nodi + correnti sorgenti di tensione)
        G: Matrice di conduttanza
        B: Matrice di accoppiamento per sorgenti di tensione
        C: Matrice di accoppiamento per correnti di sorgenti di tensione
        D: Matrice nulla per sorgenti di tensione

        Per l'analisi transitoria, i componenti dinamici (C, L) sono modellati come
        resistenze equivalenti e sorgenti di corrente/tensione equivalenti.
        Usiamo l'integrazione trapezoidale implicita.
        """
        
        # Inizializza le matrici e il vettore b
        # A_matrix = np.zeros((self._total_equations, self._total_equations))
        # b_vector = np.zeros(self._total_equations)
        
        # G matrix (conduttanze) e b vector (correnti di sorgente, ecc.)
        G_current = np.zeros((self._num_nodes, self._num_nodes))
        b_current = np.zeros(self._num_nodes)
        
        # B matrix (voltage sources to node equations)
        B_current = np.zeros((self._num_nodes, self._num_voltage_sources))
        
        # C matrix (currents through voltage sources to current equations)
        C_current = np.zeros((self._num_voltage_sources, self._num_nodes))
        
        # D matrix (voltage sources to voltage source equations - typically zero)
        D_current = np.zeros((self._num_voltage_sources, self._num_voltage_sources))
        
        # V_s_vector (for sources in voltage source equations)
        V_s_vector = np.zeros(self._num_voltage_sources)


        # Assembliamo le matrici per i componenti lineari e non dinamici
        for comp in self.circuit.components:
            n1_idx = self._get_node_index(comp.node1) if hasattr(comp, 'node1') else -1
            n2_idx = self._get_node_index(comp.node2) if hasattr(comp, 'node2') else -1

            if isinstance(comp, Resistor):
                g = 1.0 / comp.resistance
                if n1_idx != -1: G_current[n1_idx, n1_idx] += g
                if n2_idx != -1: G_current[n2_idx, n2_idx] += g
                if n1_idx != -1 and n2_idx != -1:
                    G_current[n1_idx, n2_idx] -= g
                    G_current[n2_idx, n1_idx] -= g

            elif isinstance(comp, VoltageSource):
                vs_idx = self._get_voltage_source_index(comp.name) + self._num_nodes
                
                # Aggiungi termini per la legge di Kirchhoff delle correnti al nodo
                if n1_idx != -1:
                    B_current[n1_idx, vs_idx - self._num_nodes] += 1
                    C_current[vs_idx - self._num_nodes, n1_idx] += 1 # Corrente che esce dal nodo 1

                if n2_idx != -1:
                    B_current[n2_idx, vs_idx - self._num_nodes] -= 1
                    C_current[vs_idx - self._num_nodes, n2_idx] -= 1 # Corrente che entra nel nodo 2

                # Equazione della sorgente di tensione V_n1 - V_n2 = V_source
                V_s_vector[vs_idx - self._num_nodes] = comp.get_voltage(time)

            elif isinstance(comp, CurrentSource):
                # Aggiungi la corrente al vettore b
                if n1_idx != -1: b_current[n1_idx] -= comp.current # Corrente che esce dal nodo n1
                if n2_idx != -1: b_current[n2_idx] += comp.current # Corrente che entra nel nodo n2

            elif is_transient and isinstance(comp, Capacitor):
                # Modello per l'integrazione trapezoidale implicita di un condensatore
                # I_C(t) = (2C/dt) * (V_C(t) - V_C(t-dt)) - I_C(t-dt)
                # O, in termini MNA: I_C = G_eq * V_C - I_eq
                # G_eq = 2C/dt
                # I_eq = (2C/dt)*V_C(t-dt) + I_C(t-dt)
                
                g_eq = 2.0 * comp.capacitance / dt
                
                # V_C(t) = V_n1(t) - V_n2(t)
                if n1_idx != -1: G_current[n1_idx, n1_idx] += g_eq
                if n2_idx != -1: G_current[n2_idx, n2_idx] += g_eq
                if n1_idx != -1 and n2_idx != -1:
                    G_current[n1_idx, n2_idx] -= g_eq
                    G_current[n2_idx, n1_idx] -= g_eq

                # I_eq parte
                # Calcola V_C(t-dt) = V_n1(t-dt) - V_n2(t-dt)
                V_C_prev = 0.0
                if n1_idx != -1 and prev_x is not None: V_C_prev += prev_x[n1_idx]
                if n2_idx != -1 and prev_x is not None: V_C_prev -= prev_x[n2_idx]

                # I_C(t-dt) è la corrente che scorreva nel condensatore al passo precedente.
                # Questa corrente era l'ultima corrente calcolata per questo condensatore.
                # Per la prima iterazione o se non è tracciata: si può approssimare con 0 o da IC.
                # Per implementazioni più robuste, i componenti dinamici tracciano il loro stato interno.
                
                # Qui usiamo un'approssimazione che la corrente al passo precedente è data da
                # (prev_x[n1_idx] - prev_x[n2_idx]) * g_eq_prev (se avessimo un g_eq_prev)
                # Un approccio comune è usare la corrente calcolata al passo precedente.
                # Per ora, semplifichiamo I_C(t-dt) a 0 per il primo passo, o dal valore di prev_x.
                
                # Aggiungi I_eq al vettore b
                # I_eq = g_eq * V_C(t-dt) + I_C(t-dt) (se tracciata)
                # Semplifichiamo a I_eq = g_eq * V_C(t-dt) per iniziare
                i_eq = g_eq * V_C_prev
                
                if n1_idx != -1: b_current[n1_idx] -= i_eq # Corrente che esce da n1
                if n2_idx != -1: b_current[n2_idx] += i_eq # Corrente che entra in n2

            # elif is_transient and isinstance(comp, Inductor):
            #     # Simile al condensatore, ma con modello I_L = V_L / R_eq + I_L_eq
            #     # dove R_eq = dt / (2L)
            #     # V_L(t) = V_n1(t) - V_n2(t)
            #     # I_L_eq = (V_L(t-dt) / R_eq) + I_L(t-dt)
            #     r_eq = dt / (2.0 * comp.inductance)
            #     g_eq = 1.0 / r_eq
                
            #     # V_L(t) = V_n1(t) - V_n2(t)
            #     if n1_idx != -1: G_current[n1_idx, n1_idx] += g_eq
            #     if n2_idx != -1: G_current[n2_idx, n2_idx] += g_eq
            #     if n1_idx != -1 and n2_idx != -1:
            #         G_current[n1_idx, n2_idx] -= g_eq
            #         G_current[n2_idx, n1_idx] -= g_eq

            #     # I_L_eq parte
            #     V_L_prev = 0.0
            #     if n1_idx != -1 and prev_x is not None: V_L_prev += prev_x[n1_idx]
            #     if n2_idx != -1 and prev_x is not None: V_L_prev -= prev_x[n2_idx]
                
            #     # Assumiamo che la corrente precedente sia 0 o calcolata al passo precedente
            #     i_L_prev = 0 # comp.get_previous_current() se fosse tracciata nel componente
            #     i_eq = g_eq * V_L_prev + i_L_prev
                
            #     if n1_idx != -1: b_current[n1_idx] -= i_eq
            #     if n2_idx != -1: b_current[n2_idx] += i_eq
            
            elif isinstance(comp, OpAmp):
                # Modellazione di un Op-Amp ideale o quasi ideale (V_diff = 0, I_in = 0)
                # Un Op-Amp ha 3 nodi (pin invertente, non-invertente, output)
                in_inv_idx = self._get_node_index(comp.nodes['in_inv'])
                in_non_inv_idx = self._get_node_index(comp.nodes['in_non_inv'])
                out_idx = self._get_node_index(comp.nodes['out'])
                
                # Aggiungi una riga per la tensione differenziale (V+ - V- = 0)
                # Questo richiede una corrente ausiliaria attraverso l'Op-Amp o un modello semplificato
                # Il modo più comune per gli Op-Amp ideali è aggiungere un'equazione extra (come per Vs)
                # che impone V_in_non_inv - V_in_inv = 0
                
                # Creiamo una "sorgente di tensione ausiliaria" per l'OpAmp ideale.
                # Ogni OpAmp aggiunge una riga e una colonna alla matrice.
                # Questa riga impone V_in_non_inv - V_in_inv = 0.
                # La colonna associata è la corrente che scorre nel pin d'uscita dell'OpAmp.
                
                opamp_aux_idx = self._num_nodes + self._num_voltage_sources + self.circuit.nonlinear_components.index(comp)

                # Equazione V_in_non_inv - V_in_inv = 0
                if in_non_inv_idx != -1: A_matrix[opamp_aux_idx, in_non_inv_idx] += 1
                if in_inv_idx != -1: A_matrix[opamp_aux_idx, in_inv_idx] -= 1
                
                # La corrente di output dell'OpAmp influisce sul nodo di output dell'OpAmp
                # A_matrix[out_idx, opamp_aux_idx] += 1 (questo è per il contributo della corrente)
                # In un OpAmp ideale, l'impedenza di uscita è 0. Quindi la corrente di output non è limitata.
                # L'OpAmp dovrebbe essere risolto come un elemento "extra" di Vs (con una corrente ausiliaria)
                
                # Questo modello ideale necessita di espandere le matrici G e B/C/D dinamicamente,
                # o trattare gli OpAmp come sorgenti controllate di tensione (VCVS con guadagno molto alto).
                
                # L'approccio più robusto per Op-Amp è modellarli come VCVS (Voltage Controlled Voltage Source)
                # con un guadagno molto alto, e poi usare la risoluzione non lineare per OpAmp reali.
                # Per Op-Amp ideali, l'approccio è aggiungere una Vs_aux e una I_Vs_aux per ogni OpAmp,
                # ma la Vs_aux è controllata da V_diff = 0.
                
                # Per ora, se OpAmp viene usato, assumeremo che il suo comportamento ideale (Virtual Short)
                # sia gestito implicitamente dal circuito (es. con feedback esterno).
                # L'OpAmp viene trattato come un componente non lineare per fsolve.
                pass # La logica dell'OpAmp verrà gestita nella funzione di sistema non lineare

        # Crea la matrice MNA completa
        A_matrix = np.block([
            [G_current, B_current],
            [C_current, D_current]
        ])
        
        b_vector = np.concatenate((b_current, V_s_vector))

        return A_matrix, b_vector

    def _system_equations(self, x, A_linear, b_linear, dt, prev_x, time):
        """
        Sistema di equazioni MNA, inclusi i componenti non lineari.
        x: vettore delle incognite (tensioni ai nodi + correnti sorgenti di tensione)
        A_linear, b_linear: matrici MNA per i componenti lineari
        """
        
        # Inizializza il vettore delle equazioni residuali F(x)
        F = np.dot(A_linear, x) - b_linear
        
        # Aggiungi i contributi delle non linearità al vettore F
        # Ogni componente non lineare contribuisce con correnti o tensioni.
        # Per fsolve, vogliamo F(x) = 0.
        
        for comp in self.circuit.nonlinear_components:
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
                # La corrente di griglia è assunta zero in un triodo ideale, ma potrebbe essere modellata.
            
            elif isinstance(comp, OpAmp):
                # La modellazione dell'Op-Amp come non lineare per fsolve.
                # Assumiamo un comportamento ideale con guadagno molto alto, che porta a V_diff = 0.
                # E current_input = 0.
                in_inv_idx = self._get_node_index(comp.nodes['in_inv'])
                in_non_inv_idx = self._get_node_index(comp.nodes['in_non_inv'])
                out_idx = self._get_node_index(comp.nodes['out'])
                
                V_in_inv = x[in_inv_idx] if in_inv_idx != -1 else 0.0
                V_in_non_inv = x[in_non_inv_idx] if in_non_inv_idx != -1 else 0.0
                V_out = x[out_idx] if out_idx != -1 else 0.0

                # L'equazione che l'OpAmp impone (Virtual Short per feedback negativo)
                # Se il circuito non ha feedback negativo, l'OpAmp satura.
                # Una soluzione è modellare con un guadagno molto alto: V_out = A * (V_non_inv - V_inv)
                # E poi aggiungere limitazioni di saturazione.
                
                # Per fsolve, dobbiamo esprimere le correnti che l'OpAmp inietta/assorbe dai nodi.
                # Per un OpAmp ideale con feedback negativo: V_in_inv = V_in_non_inv
                # Per OpAmp reali con saturazione:
                V_diff = V_in_non_inv - V_in_inv
                
                # Aggiungiamo un termine al vettore F che costringe V_out a essere in relazione con V_diff
                # Usiamo una grande conduttanza (o corrente) per imporre V_diff = 0
                # Questo è un modo per "simulare" il virtual short per l'OpAmp ideale.
                # Il problema con un guadagno infinito è che fsolve può avere difficoltà.
                # È meglio che la _build_linear_system imponga V_in_inv = V_in_non_inv come una Vs ausiliaria.
                # Per il momento, se l'OpAmp non ha una corrente di output da aggiungere a F, non fa nulla qui.
                pass # L'OpAmp dovrebbe essere gestito come un'equazione extra nella matrice A.

        return F

    def solve_dc(self, initial_guess=None):
        """
        Risolve il circuito in DC (stato stazionario).
        Per i componenti dinamici (C, L): C = circuito aperto, L = corto circuito.
        """
        # Preparo le matrici MNA per lo stato DC (Capacitori aperti, Induttori corti)
        # Assicurati che i componenti C e L siano ignorati o trattati come aperti/corti
        # nel metodo _build_linear_system quando is_transient=False.
        
        # Inizializza x_guess
        if initial_guess is None:
            initial_guess = np.zeros(self._total_equations + self._num_nonlinear_components_extra_eq)
        else:
            # Assicurati che initial_guess abbia la dimensione corretta
            if len(initial_guess) != (self._total_equations + self._num_nonlinear_components_extra_eq):
                 print(f"Warning: Initial guess size mismatch. Expected {self._total_equations + self._num_nonlinear_components_extra_eq}, got {len(initial_guess)}. Resetting to zeros.")
                 initial_guess = np.zeros(self._total_equations + self._num_nonlinear_components_extra_eq)

        # Build linear part of MNA (independent of non-linearities)
        A_linear, b_linear = self._build_linear_system(initial_guess, None, None, None, None, is_transient=False)
        
        try:
            # Wrap _system_equations to pass fixed arguments
            def func(x):
                return self._system_equations(x, A_linear, b_linear, dt=0, prev_x=None, time=0)

            # Solve the nonlinear system using fsolve
            solution, info, ier, msg = fsolve(func, initial_guess, full_output=True)

            if ier == 1:
                # Extract node voltages and voltage source currents
                node_voltages = solution[:self._num_nodes]
                source_currents = solution[self._num_nodes:self._num_nodes + self._num_voltage_sources]
                return node_voltages, source_currents
            else:
                print(f"DC Convergence failed: {msg}")
                return None, None
        except Exception as e:
            print(f"Error during DC simulation: {e}")
            return None, None

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
        
        # Per ora, aggiungiamo lo spazio per le equazioni degli OpAmp qui.
        # Questa è una stima, deve essere calcolata correttamente.
        self._num_nonlinear_components_extra_eq = len(self.circuit.nonlinear_components) # Assunzione per OpAmp, Diodi, Mosfet, Triodi
        
        # Inizializzazione della soluzione al tempo t-dt
        # Prima, calcola il punto di riposo DC come condizione iniziale.
        prev_x_solution, _ = self.solve_dc()
        if prev_x_solution is None:
            print("Failed to converge DC solution for initial conditions. Aborting transient.")
            return None, None

        # Aggiungi correnti iniziali delle sorgenti di tensione e degli elementi non lineari
        # Se prev_x_solution è solo V_nodes, dobbiamo espanderlo per le correnti delle sorgenti di tensione
        # e per le variabili ausiliarie degli OpAmp (se modelleli così).
        
        # Inizializza prev_x (soluzione completa al tempo t-dt)
        # La soluzione DC fornisce solo tensioni ai nodi e correnti delle sorgenti di tensione.
        # Dobbiamo inizializzare le variabili extra per le non linearità se OpAmp le usa.
        initial_solution_size = self._num_nodes + self._num_voltage_sources + self._num_nonlinear_components_extra_eq
        
        # Pre-fill with DC voltages and currents
        x_prev = np.zeros(initial_solution_size)
        x_prev[:self._num_nodes] = prev_x_solution[:self._num_nodes]
        # Se prev_x_solution contiene anche le correnti delle Vs, le copiamo.
        if len(prev_x_solution) > self._num_nodes:
            x_prev[self._num_nodes : self._num_nodes + self._num_voltage_sources] = prev_x_solution[self._num_nodes : self._num_nodes + self._num_voltage_sources]

        # Vettore per salvare i risultati
        results = {node_name: [] for node_name in (output_nodes if output_nodes else self._node_map.keys())}
        results['time'] = []
        
        # Assicurati che la DelayLine sia resettata prima della simulazione
        if self.delay_line_instance:
            self.delay_line_instance.reset()
            # Imposta il tempo di ritardo iniziale, se applicabile al circuito
            # self.delay_line_instance.set_delay_time(self.circuit.delay_time_seconds) # Esempio, dipende da come lo connetti


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
                delay_input_voltage = x_prev[self._get_node_index('delay_input')] # Usa la tensione del passo precedente
                delayed_voltage = self.delay_line_instance.update(delay_input_voltage, dt)
                
                # Ora, aggiorna la sorgente di tensione "virtuale" che rappresenta l'uscita del delay
                # Il nome della sorgente deve corrispondere a quella nel circuito
                delay_output_source_comp = next((c for c in self.circuit.components if c.name == "V_delay_output_source"), None)
                if delay_output_source_comp:
                    delay_output_source_comp.set_voltage(delayed_voltage) # Imposta la tensione della sorgente

            # Build linear part of MNA (independent of non-linearities) for this timestep
            A_linear, b_linear = self._build_linear_system(x_prev, None, None, None, None, is_transient=True, dt=dt, prev_x=x_prev, time=t)

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
