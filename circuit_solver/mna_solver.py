import numpy as np
from circuit_solver.circuit import Circuit
from components.resistor import Resistor
from components.capacitor import Capacitor
from components.inductor import Inductor
from components.voltage_source import VoltageSource
from components.current_source import CurrentSource
from components.diode import Diode

class MnaSolver:
    def __init__(self, circuit: Circuit):
        self.circuit = circuit
        self.num_total_equations = circuit.get_num_total_equations()
        
        # Variabili per l'integrazione del dominio del tempo
        self.prev_solution = np.zeros(self.num_total_equations)
        self.initial_solution_guess = np.zeros(self.num_total_equations)

        # Memorizza lo stato precedente per i componenti dinamici
        # (v_C_prev, i_L_prev). Questo è fondamentale per il metodo trapezoidale.
        self.dynamic_component_states = {}
        for comp in circuit.components:
            if isinstance(comp, (Capacitor, Inductor)):
                # Inizializza con zero o valori appropriati se noti
                self.dynamic_component_states[comp.name] = {
                    'v_prev': 0.0,
                    'i_prev': 0.0
                }
            # Se hai SpeakerDriver o RibbonTweeter, inizializza anche il loro stato meccanico
            # elif isinstance(comp, (SpeakerDriver, RibbonTweeter)):
            #     self.dynamic_component_states[comp.name] = {
            #         'v_Le_prev': 0.0, 'i_Le_prev': 0.0,
            #         'v_Lmech_prev': 0.0, 'i_Lmech_prev': 0.0,
            #         'v_Cmech_prev': 0.0, 'i_Cmech_prev': 0.0
            #     }


    def solve_dc(self, max_iterations: int = 100, tolerance: float = 1e-6):
        """
        Risolve il circuito per lo stato DC (regime stazionario).
        Per i componenti dinamici (C, L), C diventa un aperto e L un corto.
        Per i componenti non lineari, usa Newton-Raphson.
        """
        print("Starting DC analysis...")
        current_solution_guess = np.zeros(self.num_total_equations) # Inizia da zero
        
        for iteration in range(max_iterations):
            A_matrix = np.zeros((self.num_total_equations, self.num_total_equations))
            B_vector = np.zeros(self.num_total_equations)

            # Contributi dei componenti lineari (R) e sorgenti
            for component in self.circuit.components:
                if isinstance(component, Resistor):
                    n1, n2 = component.node_ids
                    G = 1.0 / component.resistance
                    if n1 != 0: A_matrix[n1, n1] += G
                    if n2 != 0: A_matrix[n2, n2] += G
                    if n1 != 0 and n2 != 0:
                        A_matrix[n1, n2] -= G
                        A_matrix[n2, n1] -= G
                elif isinstance(component, VoltageSource):
                    n_plus, n_minus = component.node_ids
                    vs_idx = component.current_index
                    
                    if n_plus != 0: A_matrix[n_plus, vs_idx] += 1
                    if n_minus != 0: A_matrix[n_minus, vs_idx] -= 1
                    
                    if n_plus != 0: A_matrix[vs_idx, n_plus] += 1
                    if n_minus != 0: A_matrix[vs_idx, n_minus] -= 1
                    B_vector[vs_idx] = component.get_voltage(0) # DC, time = 0
                elif isinstance(component, CurrentSource):
                    n_plus, n_minus = component.node_ids
                    current_val = component.get_current(0) # DC, time = 0
                    if n_plus != 0: B_vector[n_plus] -= current_val
                    if n_minus != 0: B_vector[n_minus] += current_val
                # I condensatori sono aperti in DC, gli induttori sono corti in DC (sono già gestiti dalla matrice A_matrix per i resistori)
                # La gestione dei componenti dinamici in DC è una semplificazione del metodo trapezoidale con dt -> inf

            # Contributi dei componenti non lineari (Newton-Raphson)
            # F(x) = J * dx
            F_vector = np.zeros(self.num_total_equations) # Vettore residuo
            J_matrix = np.zeros((self.num_total_equations, self.num_total_equations)) # Matrice Jacobiana

            for component in self.circuit.components:
                if isinstance(component, Diode):
                    anode_id, cathode_id = component.node_ids
                    Vd = current_solution_guess[anode_id] - current_solution_guess[cathode_id]
                    
                    # Contributo al vettore F (corrente di Shockley)
                    Id = component.calculate_current(Vd)
                    if anode_id != 0: F_vector[anode_id] -= Id
                    if cathode_id != 0: F_vector[cathode_id] += Id

                    # Contributo alla matrice Jacobiana (conduttanza differenziale)
                    Gd = component.calculate_conductance(Vd)
                    if anode_id != 0: J_matrix[anode_id, anode_id] += Gd
                    if cathode_id != 0: J_matrix[cathode_id, cathode_id] += Gd
                    if anode_id != 0 and cathode_id != 0:
                        J_matrix[anode_id, cathode_id] -= Gd
                        J_matrix[cathode_id, anode_id] -= Gd
                # Aggiungi qui altri componenti non lineari (Triode, MOSFET, ecc.)
                # e i loro contributi a F_vector e J_matrix.

            # Assembla il sistema per Newton-Raphson: J * delta_x = -F_linear - F_nonlinear
            # Dove F_linear è B_vector (il lato RHS dai componenti lineari) e F_nonlinear è F_vector (dai componenti non lineari).
            # La matrice MNA finale per Newton-Raphson è MNA_linear + J_nonlinear
            # Il vettore RHS finale per Newton-Raphson è MNA_linear * current_solution_guess - F_nonlinear
            
            # Matrice MNA totale per il passo di Newton-Raphson
            # A_total = A_matrix + J_matrix
            
            # Vettore RHS per Newton-Raphson: F(x_k) - J(x_k) * x_k
            # Per una formulazione standard F(x) = 0: J_k * delta_x_k = -F(x_k)
            # Dove F(x_k) è la somma delle correnti che lasciano i nodi
            
            # Costruiamo il vettore F(x) residuo totale (correnti ai nodi = 0)
            residual_vector = np.zeros(self.num_total_equations)
            
            # Aggiungi i contributi lineari (KCL per R, L, C nel dominio del tempo, e sorgenti)
            # La logica del get_stamps è per l'assemblaggio di G e B (stamps)
            # In Newton-Raphson, ricalcoliamo il residuo e la Jacobiana
            
            # Per i componenti lineari (Resistor, L/C in DC/TRAP)
            # R: G * (Vn1 - Vn2)
            # Vs: V_diff - V_source = 0; I_vs
            # Is: I_source
            
            # Ricostruiamo il residuo e la Jacobiana in ogni iterazione
            G_linear = np.zeros((self.num_total_equations, self.num_total_equations))
            b_linear = np.zeros(self.num_total_equations)
            
            for component in self.circuit.components:
                # Per il regime DC, i C sono aperti (G=0, I=0) e gli L sono corti (G=inf, ma gestito implicitamente)
                # Usiamo solo i contributi DC dei componenti dinamici per semplicità in DC solve.
                if isinstance(component, Resistor):
                    n1, n2 = component.node_ids
                    G = 1.0 / component.resistance
                    if n1 != 0: G_linear[n1, n1] += G
                    if n2 != 0: G_linear[n2, n2] += G
                    if n1 != 0 and n2 != 0:
                        G_linear[n1, n2] -= G
                        G_linear[n2, n1] -= G
                elif isinstance(component, VoltageSource):
                    n_plus, n_minus = component.node_ids
                    vs_idx = component.current_index
                    if n_plus != 0: G_linear[n_plus, vs_idx] += 1
                    if n_minus != 0: G_linear[n_minus, vs_idx] -= 1
                    if n_plus != 0: G_linear[vs_idx, n_plus] += 1
                    if n_minus != 0: G_linear[vs_idx, n_minus] -= 1
                    b_linear[vs_idx] = component.get_voltage(0)
                elif isinstance(component, CurrentSource):
                    n_plus, n_minus = component.node_ids
                    current_val = component.get_current(0)
                    if n_plus != 0: b_linear[n_plus] -= current_val
                    if n_minus != 0: b_linear[n_minus] += current_val

            # Calcola il vettore residuo F(x_k)
            residual_vector = G_linear @ current_solution_guess - b_linear # KCL equations, and voltage source equations

            # Aggiungi i contributi non lineari al residuo e alla Jacobiana
            J_nonlinear = np.zeros((self.num_total_equations, self.num_total_equations))
            for component in self.circuit.components:
                if isinstance(component, Diode):
                    anode_id, cathode_id = component.node_ids
                    Vd = current_solution_guess[anode_id] - current_solution_guess[cathode_id]
                    
                    Id = component.calculate_current(Vd)
                    Gd = component.calculate_conductance(Vd)

                    # Aggiungi la corrente del diodo al residuo
                    if anode_id != 0: residual_vector[anode_id] += Id
                    if cathode_id != 0: residual_vector[cathode_id] -= Id

                    # Aggiungi la conduttanza del diodo alla Jacobiana
                    if anode_id != 0: J_nonlinear[anode_id, anode_id] += Gd
                    if cathode_id != 0: J_nonlinear[cathode_id, cathode_id] += Gd
                    if anode_id != 0 and cathode_id != 0:
                        J_nonlinear[anode_id, cathode_id] -= Gd
                        J_nonlinear[cathode_id, anode_id] -= Gd
                # Aggiungi qui la logica per gli altri componenti non lineari.
                # Per Triode, MOSFET, ecc., la logica sarebbe simile, ma con più terminali e derivate incrociate.

            # Matrice Jacobiana totale per Newton-Raphson: J_total = G_linear + J_nonlinear
            J_total = G_linear + J_nonlinear

            # Risolvi il sistema per il passo di aggiornamento
            # J_total * delta_x = -residual_vector
            try:
                delta_x = np.linalg.solve(J_total, -residual_vector)
            except np.linalg.LinAlgError:
                print("Singular matrix encountered during DC solve. Check circuit topology or component values.")
                return None

            # Aggiorna la soluzione corrente
            current_solution_guess += delta_x

            # Controlla la convergenza
            if np.linalg.norm(delta_x) < tolerance:
                print(f"DC converged in {iteration + 1} iterations.")
                self.prev_solution = current_solution_guess
                return current_solution_guess
        
        print(f"DC did not converge after {max_iterations} iterations.")
        return None


    def simulate_transient(self, start_time: float, end_time: float, time_step: float,
                           max_iterations: int = 100, tolerance: float = 1e-6):
        """
        Esegue una simulazione transitoria del circuito.
        Args:
            start_time (float): Tempo di inizio della simulazione.
            end_time (float): Tempo di fine della simulazione.
            time_step (float): Passo temporale per l'integrazione.
        Returns:
            tuple: (times, solution_history)
        """
        print("Starting transient analysis...")
        times = np.arange(start_time, end_time + time_step, time_step)
        solution_history = []

        # Risolvi lo stato iniziale DC se non già fatto
        if np.all(self.prev_solution == 0):
            self.prev_solution = self.solve_dc()
            if self.prev_solution is None:
                print("DC solve failed, cannot proceed with transient simulation.")
                return [], []
        
        # Inizializza gli stati dinamici con i valori DC (o 0)
        # Per C, tensione iniziale è la tensione nodale.
        # Per L, corrente iniziale è la corrente Vs_idx
        for comp in self.circuit.components:
            if isinstance(comp, Capacitor):
                n1, n2 = comp.node_ids
                self.dynamic_component_states[comp.name]['v_prev'] = self.prev_solution[n1] - self.prev_solution[n2]
                self.dynamic_component_states[comp.name]['i_prev'] = 0.0 # Corrente iniziale non conosciuta direttamente
            elif isinstance(comp, Inductor):
                n1, n2 = comp.node_ids
                # La corrente dell'induttore in DC è la corrente attraverso il corto, difficile da estrarre
                # Potresti aver bisogno di una variabile di corrente esplicita per l'induttore in MNA per questo
                # Per semplicità qui, assumiamo 0 o un valore calcolato se L è un corto in DC.
                self.dynamic_component_states[comp.name]['v_prev'] = 0.0
                self.dynamic_component_states[comp.name]['i_prev'] = 0.0


        for t_idx, current_time in enumerate(times):
            print(f"Time step: {current_time:.6f}s")
            
            # Inizializza la guess per il passo temporale corrente con la soluzione precedente
            current_solution_guess = np.copy(self.prev_solution)
            
            for iteration in range(max_iterations):
                # Matrice Jacobiana totale e vettore residuo per Newton-Raphson
                J_total = np.zeros((self.num_total_equations, self.num_total_equations))
                residual_vector = np.zeros(self.num_total_equations)

                # Contributi dei componenti al sistema
                for component in self.circuit.components:
                    # Contributi dei componenti lineari (Resistori)
                    if isinstance(component, Resistor):
                        n1, n2 = component.node_ids
                        G = 1.0 / component.resistance
                        if n1 != 0: J_total[n1, n1] += G
                        if n2 != 0: J_total[n2, n2] += G
                        if n1 != 0 and n2 != 0:
                            J_total[n1, n2] -= G
                            J_total[n2, n1] -= G
                        
                        # Contributo al residuo: G * (Vn1 - Vn2)
                        if n1 != 0: residual_vector[n1] += G * (current_solution_guess[n1] - current_solution_guess[n2])
                        if n2 != 0: residual_vector[n2] += G * (current_solution_guess[n2] - current_solution_guess[n1])

                    # Contributi dei componenti dinamici (Capacitori, Induttori) - Metodo Trapezoidale
                    elif isinstance(component, Capacitor):
                        n1, n2 = component.node_ids
                        C = component.capacitance
                        v_C_prev = self.dynamic_component_states[component.name]['v_prev']
                        i_C_prev = self.dynamic_component_states[component.name]['i_prev']

                        # Conduttanza equivalente: G_eq = 2*C/dt
                        G_eq = 2.0 * C / time_step
                        # Corrente equivalente: i_eq = G_eq * V_C_prev + i_C_prev
                        i_eq = G_eq * v_C_prev + i_C_prev

                        if n1 != 0: J_total[n1, n1] += G_eq
                        if n2 != 0: J_total[n2, n2] += G_eq
                        if n1 != 0 and n2 != 0:
                            J_total[n1, n2] -= G_eq
                            J_total[n2, n1] -= G_eq
                        
                        # Aggiungi corrente equivalente al residuo (come sorgente)
                        # Current for Trapezoidal: i_c = G_eq * V_c(curr) - i_eq
                        # Residual at node 1: ... + i_c
                        # Residual at node 2: ... - i_c
                        residual_term_from_capacitor = G_eq * (current_solution_guess[n1] - current_solution_guess[n2]) - i_eq
                        if n1 != 0: residual_vector[n1] += residual_term_from_capacitor
                        if n2 != 0: residual_vector[n2] -= residual_term_from_capacitor

                    elif isinstance(component, Inductor):
                        n1, n2 = component.node_ids
                        L = component.inductance
                        v_L_prev = self.dynamic_component_states[component.name]['v_prev']
                        i_L_prev = self.dynamic_component_states[component.name]['i_prev']

                        # Conduttanza equivalente: G_eq = dt / (2*L)
                        G_eq = time_step / (2.0 * L)
                        # Tensione equivalente: v_eq = i_L_prev * (2*L/dt) + v_L_prev
                        v_eq = i_L_prev * (2.0 * L / time_step) + v_L_prev

                        if n1 != 0: J_total[n1, n1] += G_eq
                        if n2 != 0: J_total[n2, n2] += G_eq
                        if n1 != 0 and n2 != 0:
                            J_total[n1, n2] -= G_eq
                            J_total[n2, n1] -= G_eq

                        # Corrente equivalente da aggiungere al residuo
                        # i_L = G_eq * (V_curr) - G_eq * V_eq
                        residual_term_from_inductor = G_eq * (current_solution_guess[n1] - current_solution_guess[n2]) - G_eq * v_eq
                        if n1 != 0: residual_vector[n1] += residual_term_from_inductor
                        if n2 != 0: residual_vector[n2] -= residual_term_from_inductor


                    # Contributi delle sorgenti (indipendenti)
                    elif isinstance(component, VoltageSource):
                        n_plus, n_minus = component.node_ids
                        vs_idx = component.current_index

                        # Equazione KCL per la corrente Vs
                        if n_plus != 0: J_total[n_plus, vs_idx] += 1
                        if n_minus != 0: J_total[n_minus, vs_idx] -= 1
                        
                        # Equazione Vs: V_plus - V_minus - V_source = 0
                        if n_plus != 0: J_total[vs_idx, n_plus] += 1
                        if n_minus != 0: J_total[vs_idx, n_minus] -= 1

                        # Residuo della sorgente di tensione
                        residual_vector[vs_idx] += (current_solution_guess[n_plus] - current_solution_guess[n_minus]) - component.get_voltage(current_time)
                    
                    elif isinstance(component, CurrentSource):
                        n_plus, n_minus = component.node_ids
                        current_val = component.get_current(current_time)
                        # Contributo al residuo (corrente in uscita dal nodo)
                        if n_plus != 0: residual_vector[n_plus] += current_val
                        if n_minus != 0: residual_vector[n_minus] -= current_val

                    # Contributi dei componenti non lineari (Diodo)
                    elif isinstance(component, Diode):
                        anode_id, cathode_id = component.node_ids
                        Vd = current_solution_guess[anode_id] - current_solution_guess[cathode_id]
                        
                        Id = component.calculate_current(Vd)
                        Gd = component.calculate_conductance(Vd)

                        # Aggiungi corrente del diodo al residuo
                        if anode_id != 0: residual_vector[anode_id] += Id
                        if cathode_id != 0: residual_vector[cathode_id] -= Id

                        # Aggiungi conduttanza del diodo alla Jacobiana
                        if anode_id != 0: J_total[anode_id, anode_id] += Gd
                        if cathode_id != 0: J_total[cathode_id, cathode_id] += Gd
                        if anode_id != 0 and cathode_id != 0:
                            J_total[anode_id, cathode_id] -= Gd
                            J_total[cathode_id, anode_id] -= Gd
                    
                    # Aggiungi qui la lo
