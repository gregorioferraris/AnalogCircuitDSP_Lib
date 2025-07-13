# circuit_solver/mna_solver.py

import numpy as np
from scipy.sparse import lil_matrix # Se usi matrici sparse
from scipy.sparse.linalg import spsolve # Se usi matrici sparse
# ... importa i tuoi componenti (Resistor, VoltageSource, Triode, BJT, JFET, Diode, MOSFET, Capacitor, Inductor)
from components.resistor import Resistor
from components.voltage_source import VoltageSource
# ... altri componenti
from components.triode import Triode
from components.bjt import BJT
from components.jfet import JFET
from components.mosfet import MOSFET # Assicurati che esista e sia coerente
from components.diode import Diode

# Importa le tue funzioni Newton-Raphson e numerical_jacobian
# Assicurati che siano accessibili, magari in un modulo 'solvers' o 'utils'
# Esempio: dal file che hai appena mostrato, potresti importarle così:
from circuit_solver.nonlinear_solver import newton_raphson_solver, numerical_jacobian # <--- Assumi che tu le salvi qui

class MnaSolver:
    def __init__(self, circuit):
        self.circuit = circuit
        self.num_nodes = circuit.get_num_nodes()
        self.num_voltage_sources = circuit.get_num_voltage_sources()
        self.num_variables = self.num_nodes + self.num_voltage_sources
        self.node_map = circuit.get_node_map() # Mappa i nomi dei nodi agli indici della matrice

    def _get_node_index(self, node_name):
        return self.node_map.get(node_name, -1) # Restituisce -1 se il nodo non esiste (es. Ground)

    def solve_dc(self, initial_guess=None):
        # Prepara il guess iniziale
        if initial_guess is None:
            initial_guess = np.zeros(self.num_variables)
        elif len(initial_guess) != self.num_variables:
            raise ValueError(f"Dimensione guess iniziale ({len(initial_guess)}) non corrisponde al numero di variabili ({self.num_variables})")

        # --- Definizione delle funzioni f(x) e J(x) per Newton-Raphson ---
        def f_nonlinear_mna(x_current_guess):
            # Questa funzione costruirà il vettore F(x) = MNA_equations(x)
            # F(x) deve essere zero alla soluzione.
            
            # Matrice A (conduttanze lineari) e vettore B (sorgenti lineari)
            A_linear = lil_matrix((self.num_variables, self.num_variables), dtype=float)
            B_linear = np.zeros(self.num_variables, dtype=float)

            # Assembla le parti lineari (Resistors, Linear Sources)
            # (Codice per l'assemblaggio di A_linear e B_linear)
            for comp in self.circuit.components:
                # ... logica per resistori, sorgenti di tensione, ecc.
                # Questo è l'assemblaggio MNA tradizionale.
                if isinstance(comp, Resistor):
                    n1_idx = self._get_node_index(comp.nodes[0])
                    n2_idx = self._get_node_index(comp.nodes[1])
                    g = 1.0 / comp.resistance

                    if n1_idx != -1:
                        A_linear[n1_idx, n1_idx] += g
                    if n2_idx != -1:
                        A_linear[n2_idx, n2_idx] += g
                    if n1_idx != -1 and n2_idx != -1:
                        A_linear[n1_idx, n2_idx] -= g
                        A_linear[n2_idx, n1_idx] -= g
                elif isinstance(comp, VoltageSource):
                    # Gestione della riga extra per le sorgenti di tensione
                    pos_idx = self._get_node_index(comp.nodes['pos'])
                    neg_idx = self._get_node_index(comp.nodes['neg'])
                    source_idx = self.num_nodes + comp.index # L'indice della variabile di corrente della sorgente

                    if pos_idx != -1:
                        A_linear[pos_idx, source_idx] += 1.0
                        A_linear[source_idx, pos_idx] += 1.0
                    if neg_idx != -1:
                        A_linear[neg_idx, source_idx] -= 1.0
                        A_linear[source_idx, neg_idx] -= 1.0
                    B_linear[source_idx] += comp.voltage
                # ... altri componenti lineari come Capacitor/Inductor se in DC analysis diventano aperti/corti
                # In analisi DC, i capacitori sono circuiti aperti (nessun contributo).
                # Gli induttori sono cortocircuiti (nessun contributo alla matrice A, ma possono creare un nodo aggiuntivo o un'equazione di corrente implicita se modellati come corti).
                # Qui ci concentriamo sui componenti che hanno una presenza diretta nella matrice MNA DC.
            
            # Vettore F_nonlinear (correnti dei componenti non lineari)
            # Inizializza F_nonlinear a zero per le correnti che entrano/escono dai nodi
            F_nonlinear = np.zeros(self.num_variables, dtype=float)

            # Aggiungi i contributi dei componenti non lineari a F_nonlinear
            for comp in self.circuit.nonlinear_components:
                # Recupera le tensioni attuali dai nodi del componente dalla guess attuale
                # e calcola le correnti non lineari.
                # ATTENZIONE: Questo è il punto cruciale dove chiami i metodi dei tuoi componenti
                if isinstance(comp, Triode):
                    a_idx = self._get_node_index(comp.nodes['anode'])
                    g_idx = self._get_node_index(comp.nodes['grid'])
                    k_idx = self._get_node_index(comp.nodes['cathode'])

                    v_anode = x_current_guess[a_idx] if a_idx != -1 else 0.0
                    v_grid = x_current_guess[g_idx] if g_idx != -1 else 0.0
                    v_cathode = x_current_guess[k_idx] if k_idx != -1 else 0.0
                    
                    # Calcola Ia (corrente di placca)
                    Ia = comp.calculate_plate_current(v_grid - v_cathode, v_anode - v_cathode)

                    # Ic entra nell'anodo, esce dal catodo
                    if a_idx != -1: F_nonlinear[a_idx] += Ia
                    if k_idx != -1: F_nonlinear[k_idx] -= Ia # Corrente che esce dal catodo

                elif isinstance(comp, BJT):
                    c_idx = self._get_node_index(comp.nodes['collector'])
                    b_idx = self._get_node_index(comp.nodes['base'])
                    e_idx = self._get_node_index(comp.nodes['emitter'])

                    v_collector = x_current_guess[c_idx] if c_idx != -1 else 0.0
                    v_base = x_current_guess[b_idx] if b_idx != -1 else 0.0
                    v_emitter = x_current_guess[e_idx] if e_idx != -1 else 0.0

                    Ic_current = comp.calculate_collector_current(v_base - v_emitter, v_collector - v_emitter)
                    Ib_current = comp.calculate_base_current(Ic_current)
                    # Ie_current = Ic_current + Ib_current # Corrente che esce dall'emettitore

                    if c_idx != -1: F_nonlinear[c_idx] += Ic_current
                    if b_idx != -1: F_nonlinear[b_idx] += Ib_current
                    if e_idx != -1: F_nonlinear[e_idx] -= (Ic_current + Ib_current) # Ie esce

                elif isinstance(comp, JFET):
                    d_idx = self._get_node_index(comp.nodes['drain'])
                    g_idx = self._get_node_index(comp.nodes['gate'])
                    s_idx = self._get_node_index(comp.nodes['source'])

                    v_drain = x_current_guess[d_idx] if d_idx != -1 else 0.0
                    v_gate = x_current_guess[g_idx] if g_idx != -1 else 0.0
                    v_source = x_current_guess[s_idx] if s_idx != -1 else 0.0

                    Id_current = comp.calculate_drain_current(v_gate - v_source, v_drain - v_source)
                    
                    if d_idx != -1: F_nonlinear[d_idx] += Id_current
                    if s_idx != -1: F_nonlinear[s_idx] -= Id_current # Is esce ed è circa -Id

                elif isinstance(comp, MOSFET): # Assumiamo una classe MOSFET simile a JFET
                    d_idx = self._get_node_index(comp.nodes['drain'])
                    g_idx = self._get_node_index(comp.nodes['gate'])
                    s_idx = self._get_node_index(comp.nodes['source'])

                    v_drain = x_current_guess[d_idx] if d_idx != -1 else 0.0
                    v_gate = x_current_guess[g_idx] if g_idx != -1 else 0.0
                    v_source = x_current_guess[s_idx] if s_idx != -1 else 0.0

                    Id_current = comp.calculate_drain_current(v_gate - v_source, v_drain - v_source) # O V_ds
                    
                    if d_idx != -1: F_nonlinear[d_idx] += Id_current
                    if s_idx != -1: F_nonlinear[s_idx] -= Id_current

                elif isinstance(comp, Diode):
                    a_idx = self._get_node_index(comp.nodes['anode'])
                    c_idx = self._get_node_index(comp.nodes['cathode'])
                    
                    v_anode = x_current_guess[a_idx] if a_idx != -1 else 0.0
                    v_cathode = x_current_guess[c_idx] if c_idx != -1 else 0.0
                    
                    Id = comp.calculate_diode_current(v_anode - v_cathode)

                    if a_idx != -1: F_nonlinear[a_idx] += Id
                    if c_idx != -1: F_nonlinear[c_idx] -= Id

            # Costruisci il vettore delle equazioni f(x)
            # f(x) = A_linear @ x_current_guess + F_nonlinear - B_linear = 0
            # dove A_linear @ x_current_guess è il contributo lineare KCL/KVL
            # F_nonlinear è il contributo delle correnti non lineari
            # B_linear sono le sorgenti indipendenti
            
            # Qui A_linear è in formato lil_matrix, converti in array denso per il prodotto
            # o usa direttamente i metodi sparse se gestisci l'intera matrice A come sparsa.
            return A_linear @ x_current_guess + F_nonlinear - B_linear

        def J_nonlinear_mna(x_current_guess):
            # Questa funzione costruirà la matrice Jacobiana J(x) per f(x)
            
            # La Jacobiana complessiva è la somma della Jacobiana delle parti lineari
            # (che è semplicemente A_linear) e la Jacobiana delle parti non lineari.
            
            # Inizia con la parte lineare (A_linear è la Jacobiana dei termini lineari)
            J_mna = lil_matrix((self.num_variables, self.num_variables), dtype=float)
            
            # Assembla le parti lineari di A_linear (la stessa A_linear di f_nonlinear_mna)
            # Questo è necessario perché A_linear non è solo un termine fisso, ma una matrice di derivate
            # delle correnti lineari rispetto alle tensioni.
            for comp in self.circuit.components:
                if isinstance(comp, Resistor):
                    n1_idx = self._get_node_index(comp.nodes[0])
                    n2_idx = self._get_node_index(comp.nodes[1])
                    g = 1.0 / comp.resistance

                    if n1_idx != -1:
                        J_mna[n1_idx, n1_idx] += g
                    if n2_idx != -1:
                        J_mna[n2_idx, n2_idx] += g
                    if n1_idx != -1 and n2_idx != -1:
                        J_mna[n1_idx, n2_idx] -= g
                        J_mna[n2_idx, n1_idx] -= g
                elif isinstance(comp, VoltageSource):
                    pos_idx = self._get_node_index(comp.nodes['pos'])
                    neg_idx = self._get_node_index(comp.nodes['neg'])
                    source_idx = self.num_nodes + comp.index 
                    if pos_idx != -1:
                        J_mna[pos_idx, source_idx] += 1.0
                        J_mna[source_idx, pos_idx] += 1.0
                    if neg_idx != -1:
                        J_mna[neg_idx, source_idx] -= 1.0
                        J_mna[source_idx, neg_idx] -= 1.0
            
            # Aggiungi i contributi della Jacobiana dei componenti non lineari
            for comp in self.circuit.nonlinear_components:
                # Ottieni le tensioni dei nodi coinvolti
                component_node_order = [] # Lista dei nomi dei nodi del componente in ordine
                v_nodes_comp = [] # Lista delle tensioni assolute corrispondenti

                # Qui mappi i nomi dei nodi specifici di ciascun componente
                if isinstance(comp, Triode):
                    component_node_order = ['anode', 'grid', 'cathode']
                    v_nodes_comp = [
                        x_current_guess[self._get_node_index(comp.nodes['anode'])] if self._get_node_index(comp.nodes['anode']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['grid'])] if self._get_node_index(comp.nodes['grid']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['cathode'])] if self._get_node_index(comp.nodes['cathode']) != -1 else 0.0
                    ]
                elif isinstance(comp, BJT):
                    component_node_order = ['collector', 'base', 'emitter']
                    v_nodes_comp = [
                        x_current_guess[self._get_node_index(comp.nodes['collector'])] if self._get_node_index(comp.nodes['collector']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['base'])] if self._get_node_index(comp.nodes['base']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['emitter'])] if self._get_node_index(comp.nodes['emitter']) != -1 else 0.0
                    ]
                elif isinstance(comp, JFET):
                    component_node_order = ['drain', 'gate', 'source']
                    v_nodes_comp = [
                        x_current_guess[self._get_node_index(comp.nodes['drain'])] if self._get_node_index(comp.nodes['drain']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['gate'])] if self._get_node_index(comp.nodes['gate']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['source'])] if self._get_node_index(comp.nodes['source']) != -1 else 0.0
                    ]
                elif isinstance(comp, MOSFET):
                    component_node_order = ['drain', 'gate', 'source'] # Assumendo 3 nodi per semplicità
                    v_nodes_comp = [
                        x_current_guess[self._get_node_index(comp.nodes['drain'])] if self._get_node_index(comp.nodes['drain']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['gate'])] if self._get_node_index(comp.nodes['gate']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['source'])] if self._get_node_index(comp.nodes['source']) != -1 else 0.0
                    ]
                elif isinstance(comp, Diode):
                    component_node_order = ['anode', 'cathode']
                    v_nodes_comp = [
                        x_current_guess[self._get_node_index(comp.nodes['anode'])] if self._get_node_index(comp.nodes['anode']) != -1 else 0.0,
                        x_current_guess[self._get_node_index(comp.nodes['cathode'])] if self._get_node_index(comp.nodes['cathode']) != -1 else 0.0
                    ]
                else:
                    # Gestire altri tipi di componenti non lineari o ignorare
                    continue

                # Calcola la Jacobiana locale del componente
                J_comp = comp.calculate_jacobian_elements(*v_nodes_comp) # Passa le tensioni come argomenti separati

                # Mappa gli elementi della Jacobiana del componente nella Jacobiana globale MNA
                for row_idx_comp, row_name_comp in enumerate(component_node_order):
                    global_row_idx = self._get_node_index(comp.nodes[row_name_comp])
                    if global_row_idx == -1: # Se il nodo è Ground
                        continue 

                    for col_idx_comp, col_name_comp in enumerate(component_node_order):
                        global_col_idx = self._get_node_index(comp.nodes[col_name_comp])
                        if global_col_idx == -1: # Se il nodo è Ground
                            # Se la derivata è rispetto a Ground, non contribuisce alla matrice MNA
                            # MA le correnti che fluiscono a Ground influenzano le altre equazioni.
                            # Qui stiamo assumendo che Ground sia il nodo 0 e non abbia un'equazione KCL.
                            # Se Ground è -1, ignoriamo la colonna corrispondente alla derivata rispetto a Ground.
                            continue 
                        J_mna[global_row_idx, global_col_idx] += J_comp[row_idx_comp, col_idx_comp]
            
            return J_mna.toarray() # Restituisci la Jacobiana come array denso per newton_raphson_solver

        # --- Chiama il solutore di Newton-Raphson ---
        try:
            solution = newton_raphson_solver(f_nonlinear_mna, J_nonlinear_mna, initial_guess)
            # Estrai le tensioni ai nodi e le correnti delle sorgenti
            node_voltages = solution[:self.num_nodes]
            source_currents = solution[self.num_nodes:]
            return node_voltages, source_currents
        except RuntimeError as e:
            print(f"Errore nella risoluzione DC: {e}")
            return None, None
