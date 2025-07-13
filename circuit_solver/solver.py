# circuit_solver/solver.py

import numpy as np

# Importa le classi necessarie dal modulo circuit
from circuit_solver.circuit import Circuit

# Importa le costanti e le utility numeriche
from utils.constants import DEFAULT_SAMPLE_RATE
from utils.helpers import newton_raphson_solver, numerical_jacobian

# Importa TUTTI i tipi di componenti per poterli gestire nell'MNA
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


class CircuitSolver:
    """
    Classe responsabile della risoluzione di un circuito definito dalla classe Circuit.
    Implementa l'analisi nodale modificata (MNA) e utilizza il metodo di Newton-Raphson
    per le non linearità.
    """
    def __init__(self, circuit: Circuit, sample_rate=DEFAULT_SAMPLE_RATE):
        self.circuit = circuit
        self.sample_rate = float(sample_rate)
        self.Ts = 1.0 / self.sample_rate # Periodo di campionamento

        self.num_nodes = self.circuit.get_num_nodes() # Numero totale di nodi
        self.num_aux_vars = self._get_num_aux_vars() # Numero di variabili ausiliarie (es. per induttori)
        self.total_vars = self.num_nodes + self.num_aux_vars

        # Vettore delle incognite (tensioni nodali + correnti ausiliarie) al passo precedente
        self.x_prev = np.zeros(self.total_vars)
        self.x_prev[self.circuit.get_ground_node_id()] = 0.0 # Il nodo GND è sempre 0V

        # Dizionari per mappare i componenti che generano variabili ausiliarie ai loro indici
        self.inductor_aux_map = {} # {component_id: aux_var_index}
        # Inizializza le mappe delle variabili ausiliarie
        self._map_auxiliary_variables()

        print(f"Solutore inizializzato per circuito '{circuit.name}' con Fs={self.sample_rate}Hz.")
        print(f"Numero di nodi da risolvere: {self.num_nodes}")
        print(f"Numero di variabili ausiliarie: {self.num_aux_vars}")
        print(f"Numero totale di variabili MNA: {self.total_vars}")


    def _map_auxiliary_variables(self):
        """
        Popola i dizionari di mappatura per le variabili ausiliarie.
        Questo è fondamentale per costruire correttamente la matrice MNA.
        """
        aux_idx = 0 # Inizia dopo i nodi
        for comp in self.circuit.get_components():
            if isinstance(comp, Inductor):
                self.inductor_aux_map[comp.component_id] = self.num_nodes + aux_idx
                aux_idx += 1
            # Aggiungi altri tipi di componenti che necessitano di variabili ausiliarie qui
            # Es: sorgenti di tensione (se implementate come variabili V_aux)


    def _get_num_aux_vars(self):
        """
        Calcola il numero di variabili ausiliarie (correnti) necessarie per l'MNA.
        Es: per ogni induttore, si aggiunge una variabile per la sua corrente.
        """
        num_aux = 0
        for comp in self.circuit.get_components():
            # Ogni induttore aggiunge una variabile ausiliaria per la sua corrente.
            if isinstance(comp, Inductor):
                num_aux += 1
            # Aggiungi qui altre condizioni per componenti che introducono variabili ausiliarie
            # Es: sorgenti di tensione ideali (non in questo set di componenti)
        return num_aux


    def build_mna_system(self, x_guess, input_signal_value=0.0):
        """
        Costruisce il sistema di equazioni MNA non lineari f(x) = 0.
        Questo metodo viene passato come `func` al solutore di Newton-Raphson.

        Args:
            x_guess (np.ndarray): Il guess corrente del vettore delle incognite [V_nodes, I_aux].
            input_signal_value (float): Il valore del segnale di ingresso al campione attuale.

        Returns:
            np.ndarray: Il vettore delle funzioni f(x) = 0.
        """
        # Estrai le tensioni nodali e le correnti ausiliarie dal vettore x_guess
        node_voltages = x_guess[:self.num_nodes]
        aux_currents = x_guess[self.num_nodes:]

        # Inizializza il vettore delle equazioni f(x)
        f_vector = np.zeros(self.total_vars)

        # 1. Equazioni KCL per ogni nodo (tranne GND)
        # Queste equazioni saranno popolate dai contributi dei componenti.
        # f_vector[0] rimane 0 (riferimento GND)

        # 2. Scansiona i componenti e aggiungi i loro contributi
        for comp in self.circuit.get_components():
            n1, n2, *extra_nodes = comp.connected_nodes # Nodi a cui il componente è connesso

            # Inizializza le tensioni ai capi del componente
            v_n1 = node_voltages[n1]
            v_n2 = node_voltages[n2]
            voltage_across_comp = v_n1 - v_n2

            # Contributi per componenti passivi/non lineari a 2 terminali (Corrente(Tensione))
            if isinstance(comp, (Resistor, Diode, LED, SchottkyDiode, ZenerDiode, RectifierTube)):
                current = comp.calculate_current(voltage_across_comp)
                if n1 != self.circuit.get_ground_node_id():
                    f_vector[n1] += current
                if n2 != self.circuit.get_ground_node_id():
                    f_vector[n2] -= current
            elif isinstance(comp, LDR): # LDR è una resistenza variabile, si comporta come Resistor
                # La resistenza dell'LDR dipende dal livello di luce, che dovrebbe essere aggiornato
                # esternamente o prima di questa chiamata. Qui usiamo la resistenza corrente dell'LDR.
                current = (node_voltages[n1] - node_voltages[n2]) / comp.current_resistance # Usa la resistenza attuale dell'LDR
                if n1 != self.circuit.get_ground_node_id():
                    f_vector[n1] += current
                if n2 != self.circuit.get_ground_node_id():
                    f_vector[n2] -= current
            elif isinstance(comp, Capacitor):
                # Modellazione del condensatore con integrazione trapezoidale: I_C = C/Ts * (V_C(t) - V_C(t-Ts))
                # La corrente I_C è una funzione delle tensioni attuali e precedenti.
                v_cap_current_guess = voltage_across_comp
                v_cap_prev = comp.get_previous_voltage() # Metodo da aggiungere a Capacitor
                current_cap = comp.C / self.Ts * (v_cap_current_guess - v_cap_prev)
                if n1 != self.circuit.get_ground_node_id():
                    f_vector[n1] += current_cap
                if n2 != self.circuit.get_ground_node_id():
                    f_vector[n2] -= current_cap
            elif isinstance(comp, Inductor):
                # Modellazione dell'induttore con variabile ausiliaria (corrente dell'induttore)
                # Equazione KCL ai nodi dell'induttore: I_L entra nel nodo n1, esce dal nodo n2
                # (Assumendo che I_L sia la variabile ausiliaria)
                inductor_aux_idx = self.inductor_aux_map[comp.component_id] - self.num_nodes # Indice relativo in aux_currents
                i_L_current_guess = aux_currents[inductor_aux_idx]

                if n1 != self.circuit.get_ground_node_id():
                    f_vector[n1] += i_L_current_guess
                if n2 != self.circuit.get_ground_node_id():
                    f_vector[n2] -= i_L_current_guess

                # Equazione ausiliaria per l'induttore (KCL)
                # V_L(t) = L * (I_L(t) - I_L(t-Ts)) / Ts
                # (V_n1 - V_n2) - L/Ts * I_L(t) + L/Ts * I_L(t-Ts) = 0
                f_vector[self.num_nodes + inductor_aux_idx] = \
                    (v_n1 - v_n2) - comp.L / self.Ts * i_L_current_guess + comp.L / self.Ts * comp.get_previous_current() # Metodo da aggiungere a Inductor


            # Componenti a 3 terminali (Transistor: Drain, Gate, Source; BJT: Collector, Base, Emitter)
            elif isinstance(comp, (JFET, MOSFET)):
                n_drain, n_gate, n_source = n1, n2, extra_nodes[0] # Mappa i nodi in base alla definizione del componente
                v_ds = node_voltages[n_drain] - node_voltages[n_source]
                v_gs = node_voltages[n_gate] - node_voltages[n_source]
                id_current = comp.calculate_drain_current(v_gs, v_ds)

                if n_drain != self.circuit.get_ground_node_id():
                    f_vector[n_drain] += id_current
                if n_source != self.circuit.get_ground_node_id():
                    f_vector[n_source] -= id_current # La corrente esce dal source
                # La corrente di Gate è solitamente zero per JFET/MOSFET ideali
                # f_vector[n_gate] += I_gate (se modellata)

            elif isinstance(comp, BJT):
                n_collector, n_base, n_emitter = n1, n2, extra_nodes[0]
                v_ce = node_voltages[n_collector] - node_voltages[n_emitter]
                v_be = node_voltages[n_base] - node_voltages[n_emitter]
                ic_current = comp.calculate_collector_current(v_be, v_ce)
                ib_current = comp.calculate_base_current(v_be, v_ce) # La corrente di base non è nulla

                if n_collector != self.circuit.get_ground_node_id():
                    f_vector[n_collector] += ic_current
                if n_base != self.circuit.get_ground_node_id():
                    f_vector[n_base] += ib_current
                if n_emitter != self.circuit.get_ground_node_id():
                    f_vector[n_emitter] -= (ic_current + ib_current) # Corrente di emettitore = Ic + Ib


            # Componenti a 3 terminali (Valvole: Anode, Grid, Cathode)
            elif isinstance(comp, Triode):
                n_anode, n_grid, n_cathode = n1, n2, extra_nodes[0]
                v_a = node_voltages[n_anode] - node_voltages[n_cathode]
                v_g = node_voltages[n_grid] - node_voltages[n_cathode]
                ia_current = comp.calculate_anode_current(v_g, v_a)

                if n_anode != self.circuit.get_ground_node_id():
                    f_vector[n_anode] += ia_current
                if n_cathode != self.circuit.get_ground_node_id():
                    f_vector[n_cathode] -= ia_current
                # Corrente di Griglia (Ig) può essere modellata se Vg > 0
                # Se non modellata, Ig = 0

            # Componenti a 4/5 terminali (Pentodo: Anode, Grid1, Grid2, Cathode)
            elif isinstance(comp, Pentode):
                n_anode, n_grid1, n_grid2, n_cathode = n1, n2, extra_nodes[0], extra_nodes[1]
                v_a = node_voltages[n_anode] - node_voltages[n_cathode]
                v_g1 = node_voltages[n_grid1] - node_voltages[n_cathode]
                v_g2 = node_voltages[n_grid2] - node_voltages[n_cathode]
                ia_current = comp.calculate_anode_current(v_g1, v_g2, v_a)
                # Tipicamente, la corrente di G2 (screen grid) è una frazione di Ia, ma non inclusa in questo modello
                # Sarà necessaria una modellazione più complessa per Ig2.

                if n_anode != self.circuit.get_ground_node_id():
                    f_vector[n_anode] += ia_current
                if n_cathode != self.circuit.get_ground_node_id():
                    f_vector[n_cathode] -= ia_current


        # 3. Gestione del segnale di ingresso
        # Assumiamo che ci sia un nodo chiamato "Input" nel circuito.
        # Se c'è una sorgente di tensione ideale V_in tra "Input" e "GND",
        # l'equazione sarebbe V_Input - V_GND - V_in = 0.
        # Per implementare una sorgente di tensione, tipicamente si aggiunge una variabile ausiliaria
        # per la corrente attraverso la sorgente.
        # Per ora, una semplificazione è "forzare" la tensione sul nodo di ingresso,
        # ma questo non è MNA puro e può causare problemi.
        # Un'implementazione più robusta avrebbe una sorgente di tensione come componente.

        # Se il circuito ha un nodo 'Input', e si vuole applicare il segnale a quello
        # Potrebbe essere gestito come una sorgente di tensione controllata, ma per ora è placeholder.
        input_node_id = self.circuit.get_node_id("Input") if "Input" in self.circuit.nodes else None
        if input_node_id is not None and input_node_id != self.circuit.get_ground_node_id():
            # Questo è un trucco comune ma non rigoroso per applicare un input
            # Se vogliamo che V_input = input_signal_value, allora la sua equazione KCL sarà:
            # f_vector[input_node_id] = node_voltages[input_node_id] - input_signal_value
            # Questo trasforma l'equazione KCL di quel nodo in un'equazione V_in = V_applied.
            # La riga e colonna corrispondenti nella Jacobiana saranno [1,0,0...; 0,0,0...]
            # Questa è una forma semplificata di gestione della sorgente di tensione.
            # L'implementazione MNA completa aggiungerebbe una variabile extra per la corrente della sorgente.
            pass # L'input_signal_value sarà usato quando la MNA viene completamente riempita.

        # Nota: La riga e la colonna per il nodo di massa (GND) sono gestite nella Jacobiana.
        # L'equazione per GND è implicitamente 0 = 0 o V_GND = 0.
        f_vector[self.circuit.get_ground_node_id()] = node_voltages[self.circuit.get_ground_node_id()] # Questo assicura che V_GND = 0

        return f_vector

    def build_jacobian_matrix(self, x_guess, input_signal_value=0.0):
        """
        Costruisce la matrice Jacobiana J(x) delle equazioni MNA.
        Questo metodo viene passato come `jacobian` al solutore di Newton-Raphson.

        Args:
            x_guess (np.ndarray): Il guess corrente del vettore delle incognite [V_nodes, I_aux].
            input_signal_value (float): Il valore del segnale di ingresso al campione attuale.

        Returns:
            np.ndarray: La matrice Jacobiana J(x).
        """
        jacobian_matrix = np.zeros((self.total_vars, self.total_vars))

        # Estrai le tensioni nodali e le correnti ausiliarie dal vettore x_guess
        node_voltages = x_guess[:self.num_nodes]
        # aux_currents = x_guess[self.num_nodes:] # Non sempre necessario direttamente qui

        # 1. Contributi dalle derivate delle equazioni KCL
        for comp in self.circuit.get_components():
            n1, n2, *extra_nodes = comp.connected_nodes

            # Conduttanza (dI/dV) per componenti a 2 terminali
            # È VITALI che i tuoi componenti abbiano metodi per calcolare le loro conduttanze dinamiche!
            # Se non sono presenti, useremo una derivata numerica (meno efficiente).
            voltage_across_comp = node_voltages[n1] - node_voltages[n2]

            if isinstance(comp, Resistor):
                G = 1.0 / comp.resistance
                jacobian_matrix[n1, n1] += G
                jacobian_matrix[n1, n2] -= G
                jacobian_matrix[n2, n1] -= G
                jacobian_matrix[n2, n2] += G
            elif isinstance(comp, (Diode, LED, SchottkyDiode, ZenerDiode, RectifierTube)):
                # Questi componenti DOVREBBERO avere un metodo calculate_conductance(Vd)
                try:
                    g_d = comp.calculate_conductance(voltage_across_comp)
                except AttributeError:
                    # Fallback numerico se il metodo non esiste (per il debug, non per produzione)
                    # NOTA: Questo è MOLTO meno efficiente. Implementa calculate_conductance!
                    print(f"ATTENZIONE: {comp.__class__.__name__} non ha calculate_conductance. Usando derivata numerica.")
                    g_d = numerical_jacobian(lambda v: comp.calculate_current(v), voltage_across_comp)
                    if isinstance(g_d, np.ndarray): g_d = g_d[0] # Se numerical_jacobian restituisce un array

                jacobian_matrix[n1, n1] += g_d
                jacobian_matrix[n1, n2] -= g_d
                jacobian_matrix[n2, n1] -= g_d
                jacobian_matrix[n2, n2] += g_d
            elif isinstance(comp, LDR):
                # La LDR è una resistenza, quindi la sua conduttanza è 1/R_LDR
                G = 1.0 / comp.current_resistance # Usa la resistenza attuale dell'LDR
                jacobian_matrix[n1, n1] += G
                jacobian_matrix[n1, n2] -= G
                jacobian_matrix[n2, n1] -= G
                jacobian_matrix[n2, n2] += G
            elif isinstance(comp, Capacitor):
                # Contributo dalla discretizzazione trapezoidale: C/Ts
                G_cap = comp.C / self.Ts
                jacobian_matrix[n1, n1] += G_cap
                jacobian_matrix[n1, n2] -= G_cap
                jacobian_matrix[n2, n1] -= G_cap
                jacobian_matrix[n2, n2] += G_cap
            elif isinstance(comp, Inductor):
                # Contributi MNA dell'induttore quando la corrente è una variabile ausiliaria
                # Equazione KCL al nodo n1: ... + I_L = 0
                # Equazione KCL al nodo n2: ... - I_L = 0
                inductor_aux_mna_idx = self.inductor_aux_map[comp.component_id]

                # Derivata della KCL rispetto alla corrente ausiliaria dell'induttore
                jacobian_matrix[n1, inductor_aux_mna_idx] += 1.0
                jacobian_matrix[n2, inductor_aux_mna_idx] -= 1.0

                # Equazione ausiliaria per l'induttore: (V_n1 - V_n2) - L/Ts * I_L + L/Ts * I_L_prev = 0
                # Derivata dell'equazione ausiliaria rispetto a V_n1, V_n2 e I_L
                jacobian_matrix[inductor_aux_mna_idx, n1] += 1.0
                jacobian_matrix[inductor_aux_mna_idx, n2] -= 1.0
                jacobian_matrix[inductor_aux_mna_idx, inductor_aux_mna_idx] -= comp.L / self.Ts

            # Transistor e Valvole (contributi di transconduttanza e conduttanza)
            elif isinstance(comp, (JFET, MOSFET)):
                n_drain, n_gate, n_source = n1, n2, extra_nodes[0]
                v_ds = node_voltages[n_drain] - node_voltages[n_source]
                v_gs = node_voltages[n_gate] - node_voltages[n_source]

                # gm = d(Id)/d(Vgs), gds = d(Id)/d(Vds)
                try:
                    gm = comp.calculate_transconductance(v_gs, v_ds) # Dovrai aggiungere questo metodo
                    gds = comp.calculate_output_conductance(v_gs, v_ds) # Dovrai aggiungere questo metodo
                except AttributeError:
                    print(f"ATTENZIONE: {comp.__class__.__name__} non ha metodi di conduttanza/transconduttanza. Usando derivata numerica.")
                    # Fallback numerico più complesso per funzioni a più variabili
                    # Per ora, usiamo una semplificazione: assumiamo che sia la derivata rispetto a Vgs
                    gm = numerical_jacobian(lambda v: comp.calculate_drain_current(v, v_ds), v_gs)
                    gds = numerical_jacobian(lambda v: comp.calculate_drain_current(v_gs, v), v_ds)
                    if isinstance(gm, np.ndarray): gm = gm[0]
                    if isinstance(gds, np.ndarray): gds = gds[0]

                # Contributi all'equazione del Drain (f_Drain)
                # f_Drain = Id(Vgs, Vds)
                # d(f_Drain)/d(V_Drain) = gds
                # d(f_Drain)/d(V_Source) = -gds - gm
                # d(f_Drain)/d(V_Gate) = gm
                if n_drain != self.circuit.get_ground_node_id():
                    jacobian_matrix[n_drain, n_drain] += gds
                    jacobian_matrix[n_drain, n_source] -= (gds + gm) # Se Source è il riferimento
                    jacobian_matrix[n_drain, n_gate] += gm
                if n_source != self.circuit.get_ground_node_id():
                    # Contributi per l'equazione KCL del Source
                    # La somma delle correnti al Source = -Id. Quindi la derivata di -(Id).
                    jacobian_matrix[n_source, n_drain] -= gds
                    jacobian_matrix[n_source, n_source] += (gds + gm)
                    jacobian_matrix[n_source, n_gate] -= gm
                # Griglia (Gate) tipicamente non ha corrente in/out, quindi la riga e colonna per Gate sono 0
                # A meno che non si modelli la corrente di gate.


            elif isinstance(comp, BJT):
                n_collector, n_base, n_emitter = n1, n2, extra_nodes[0]
                v_ce = node_voltages[n_collector] - node_voltages[n_emitter]
                v_be = node_voltages[n_base] - node_voltages[n_emitter]

                # gm = d(Ic)/d(Vbe), gpi = d(Ib)/d(Vbe), ro = d(Ic)/d(Vce)
                try:
                    gm = comp.calculate_transconductance(v_be, v_ce) # dIc/dVbe
                    gpi = comp.calculate_input_conductance(v_be, v_ce) # dIb/dVbe
                    gce = comp.calculate_output_conductance(v_be, v_ce) # dIc/dVce
                except AttributeError:
                    print(f"ATTENZIONE: {comp.__class__.__name__} non ha metodi di conduttanza/transconduttanza. Usando derivata numerica.")
                    gm = numerical_jacobian(lambda v: comp.calculate_collector_current(v, v_ce), v_be)
                    gpi = numerical_jacobian(lambda v: comp.calculate_base_current(v, v_ce), v_be)
                    gce = numerical_jacobian(lambda v: comp.calculate_collector_current(v_be, v), v_ce)
                    if isinstance(gm, np.ndarray): gm = gm[0]
                    if isinstance(gpi, np.ndarray): gpi = gpi[0]
                    if isinstance(gce, np.ndarray): gce = gce[0]


                # Contributi all'equazione del Collettore (f_Collector)
                # f_Collector = Ic(Vbe, Vce)
                if n_collector != self.circuit.get_ground_node_id():
                    jacobian_matrix[n_collector, n_collector] += gce
                    jacobian_matrix[n_collector, n_emitter] -= (gm + gce) # Se Emitter è il riferimento
                    jacobian_matrix[n_collector, n_base] += gm

                # Contributi all'equazione della Base (f_Base)
                # f_Base = Ib(Vbe, Vce)
                if n_base != self.circuit.get_ground_node_id():
                    jacobian_matrix[n_base, n_base] += gpi
                    jacobian_matrix[n_base, n_emitter] -= gpi # Se Emitter è il riferimento
                    # dIb/dVce è solitamente trascurabile per BJT, ma se modellata andrebbe qui

                # Contributi all'equazione dell'Emettitore (f_Emitter)
                # f_Emitter = -(Ic + Ib)
                if n_emitter != self.circuit.get_ground_node_id():
                    jacobian_matrix[n_emitter, n_collector] -= gce
                    jacobian_matrix[n_emitter, n_base] -= gm - gpi # Attenzione ai segni!
                    jacobian_matrix[n_emitter, n_emitter] += (gm + gpi + gce) # Somma delle conduttanze


            elif isinstance(comp, Triode):
                n_anode, n_grid, n_cathode = n1, n2, extra_nodes[0]
                v_a = node_voltages[n_anode] - node_voltages[n_cathode]
                v_g = node_voltages[n_grid] - node_voltages[n_cathode]

                # gm = d(Ia)/d(Vg), rp = d(Va)/d(Ia) = 1/gp (gp = d(Ia)/d(Va))
                try:
                    gm = comp.calculate_transconductance(v_g, v_a) # dIa/dVg
                    gp = comp.calculate_plate_conductance(v_g, v_a) # dIa/dVa
                except AttributeError:
                    print(f"ATTENZIONE: {comp.__class__.__name__} non ha metodi di conduttanza/transconduttanza. Usando derivata numerica.")
                    gm = numerical_jacobian(lambda v: comp.calculate_anode_current(v, v_a), v_g)
                    gp = numerical_jacobian(lambda v: comp.calculate_anode_current(v_g, v), v_a)
                    if isinstance(gm, np.ndarray): gm = gm[0]
                    if isinstance(gp, np.ndarray): gp = gp[0]


                # Contributi all'equazione dell'Anodo (f_Anode)
                # f_Anode = Ia(Vg, Va)
                if n_anode != self.circuit.get_ground_node_id():
                    jacobian_matrix[n_anode, n_anode] += gp
                    jacobian_matrix[n_anode, n_cathode] -= (gp + gm) # Se Cathode è il riferimento
                    jacobian_matrix[n_anode, n_grid] += gm

                # Contributi all'equazione del Catodo (f_Cathode)
                # f_Cathode = -Ia(Vg, Va)
                if n_cathode != self.circuit.get_ground_node_id():
                    jacobian_matrix[n_cathode, n_anode] -= gp
                    jacobian_matrix[n_cathode, n_cathode] += (gp + gm)
                    jacobian_matrix[n_cathode, n_grid] -= gm

                # Griglia (Grid) non ha corrente (idealmente), quindi la riga e colonna per Grid sono 0
                # a meno di modellare la corrente di griglia (Ig).


            elif isinstance(comp, Pentode):
                n_anode, n_grid1, n_grid2, n_cathode = n1, n2, extra_nodes[0], extra_nodes[1]
                v_a = node_voltages[n_anode] - node_voltages[n_cathode]
                v_g1 = node_voltages[n_grid1] - node_voltages[n_cathode]
                v_g2 = node_voltages[n_grid2] - node_voltages[n_cathode]

                # gm1 = d(Ia)/d(Vg1), gm2 = d(Ia)/d(Vg2), gp = d(Ia)/d(Va)
                try:
                    gm1 = comp.calculate_transconductance_g1(v_g1, v_g2, v_a)
                    gm2 = comp.calculate_transconductance_g2(v_g1, v_g2, v_a)
                    gp = comp.calculate_plate_conductance(v_g1, v_g2, v_a)
                except AttributeError:
                    print(f"ATTENZIONE: {comp.__class__.__name__} non ha metodi di conduttanza/transconduttanza. Usando derivata numerica.")
                    gm1 = numerical_jacobian(lambda v: comp.calculate_anode_current(v, v_g2, v_a), v_g1)
                    gm2 = numerical_jacobian(lambda v: comp.calculate_anode_current(v_g1, v, v_a), v_g2)
                    gp = numerical_jacobian(lambda v: comp.calculate_anode_current(v_g1, v_g2, v), v_a)
                    if isinstance(gm1, np.ndarray): gm1 = gm1[0]
                    if isinstance(gm2, np.ndarray): gm2 = gm2[0]
                    if isinstance(gp, np.ndarray): gp = gp[0]

                # Contributi all'equazione dell'Anodo (f_Anode)
                # f_Anode = Ia(Vg1, Vg2, Va)
                if n_anode != self.circuit.get_ground_node_id():
                    jacobian_matrix[n_anode, n_anode] += gp
                    jacobian_matrix[n_anode, n_cathode] -= (gp + gm1 + gm2) # Se Cathode è il riferimento
                    jacobian_matrix[n_anode, n_grid1] += gm1
                    jacobian_matrix[n_anode, n_grid2] += gm2

                # Contributi all'equazione del Catodo (f_Cathode)
                # f_Cathode = -Ia(Vg1, Vg2, Va)
                if n_cathode != self.circuit.get_ground_node_id():
                    jacobian_matrix[n_cathode, n_anode] -= gp
                    jacobian_matrix[n_cathode, n_cathode] += (gp + gm1 + gm2)
                    jacobian_matrix[n_cathode, n_grid1] -= gm1
                    jacobian_matrix[n_cathode, n_grid2] -= gm2

                # Per Griglia1 e Griglia2, tipicamente non c'è corrente di griglia,
                # quindi le loro righe e colonne per i contributi del pentodo sono 0.
                # Se la corrente di griglia viene modellata, si aggiungono qui.

        # 2. Contributi dalle equazioni ausiliarie
        # Per gli induttori: equazione ausiliaria: (V_n1 - V_n2) - L/Ts * I_L + L/Ts * I_L_prev = 0
        # Questa equazione deriva da dV/dt = L * dI/dt => V_L = L * (I_L - I_L_prev) / Ts
        for comp in self.circuit.get_components():
            if isinstance(comp, Inductor):
                n1, n2 = comp.connected_nodes
                inductor_aux_mna_idx = self.inductor_aux_map[comp.component_id]

                # Derivata della riga aux (che è V_L - L/Ts * I_L) rispetto alle variabili
                jacobian_matrix[inductor_aux_mna_idx, n1] += 1.0
                jacobian_matrix[inductor_aux_mna_idx, n2] -= 1.0
                jacobian_matrix[inductor_aux_mna_idx, inductor_aux_mna_idx] -= comp.L / self.Ts


        # Gestione del nodo di massa (GND): Fissa V_GND a 0
        # Imposta la riga e la colonna del nodo GND (normalmente il nodo 0)
        # in modo che l'equazione sia V_GND = 0 e la derivata sia 1 rispetto a V_GND.
        gnd_id = self.circuit.get_ground_node_id()
        jacobian_matrix[gnd_id, :] = 0.0 # Cancella tutta la riga (perché V_GND è nota)
        jacobian_matrix[:, gnd_id] = 0.0 # Cancella tutta la colonna (i contributi che dipendono da V_GND sono zero)
        jacobian_matrix[gnd_id, gnd_id] = 1.0 # Imposta V_GND = 0


        # Gestione dell'input per sorgenti di tensione (se implementate come variabili ausiliarie)
        # Se c'è una sorgente di tensione V_in tra Input e GND, avrai un'equazione V_Input - V_GND = V_in.
        # Questa equazione e la sua variabile ausiliaria per la corrente della sorgente saranno gestite qui.
        input_node_id = self.circuit.get_node_id("Input") if "Input" in self.circuit.nodes else None
        if input_node_id is not None and input_node_id != gnd_id:
            # Questo è un approccio semplificato per forzare la tensione su un nodo,
            # trasformando la sua equazione KCL in una definizione di tensione.
            # Se fosse una vera sorgente di tensione ideale, la gestione sarebbe diversa.
            # Qui si modifica la Jacobiana per far sì che V_input = input_signal_value.
            # Questo approccio può essere problematico in circuiti complessi o con sorgenti di corrente.
            # La Jacobiana modificata riflette il fatto che V_input è ora una variabile nota.
            pass # Questa parte sarà implementata nel metodo solve_sample

        return jacobian_matrix


    def solve_sample(self, input_signal_value=0.0):
        """
        Risolve il circuito per un singolo campione audio.

        Args:
            input_signal_value (float): Il valore del segnale di ingresso al campione attuale.

        Returns:
            np.ndarray: Le tensioni nodali risolte per il campione attuale (solo nodi, senza aux vars).
        """
        # Il guess iniziale per Newton-Raphson è lo stato del campione precedente
        x_guess = self.x_prev.copy()

        # Applicare il segnale di ingresso.
        # Se hai un nodo "Input" designato, puoi "iniettare" l'input qui.
        # Questo è un modo comune per gestire le sorgenti di tensione in MNA:
        # Aggiungi una sorgente di tensione ideale tra il nodo di input e GND.
        # Questo aggiungerebbe una riga e una colonna alla Jacobiana, e una variabile ausiliaria
        # per la corrente attraverso la sorgente.
        # Per la nostra implementazione MNA attuale, il modo più semplice è modificare
        # il termine noto del vettore f(x) o forzare la tensione del nodo di input.

        # Metodo 1 (più robusto per MNA completa): La sorgente di tensione viene trattata come componente.
        # Metodo 2 (semplificato per i nostri placeholder): Passa input_signal_value alle funzioni build_mna_system.
        # Le funzioni build_mna_system e build_jacobian_matrix dovranno usare questo valore per modificare
        # le equazioni relative al nodo di input (es. se V_input = input_signal_value, allora f_input = V_input - input_signal_value).

        # Il metodo _solve_mna_system implementa l'algoritmo di Newton-Raphson
        try:
            # Passiamo una lambda che cattura input_signal_value
            # Questo permette a build_mna_system e build_jacobian_matrix di accedere al valore di input.
            current_x = newton_raphson_solver(
                func=lambda x_curr: self.build_mna_system(x_curr, input_signal_value),
                jacobian=lambda x_curr: self.build_jacobian_matrix(x_curr, input_signal_value),
                x0=x_guess
            )
        except RuntimeError as e:
            print(f"Errore nella risoluzione del campione: {e}")
            # In caso di non convergenza, restituisci lo stato precedente per stabilità
            return self.x_prev[:self.num_nodes] # Restituisce solo le tensioni nodali

        # Aggiorna lo stato dei componenti dinamici (condensatori, induttori, LDR)
        # e le variabili di stato del solutore per il prossimo campione.
        self._update_component_states(current_x)

        # Memorizza lo stato risolto per il prossimo passo
        self.x_prev = current_x.copy()

        # Assicurati che il nodo GND rimanga a 0V
        self.x_prev[self.circuit.get_ground_node_id()] = 0.0

        return current_x[:self.num_nodes] # Restituisci solo le tensioni nodali


    def _update_component_states(self, current_x):
        """
        Aggiorna lo stato interno dei componenti dinamici (es. C, L, LDR)
        dopo che il circuito è stato risolto per il campione attuale.
        """
        node_voltages = current_x[:self.num_nodes]
        aux_currents = current_x[self.num_nodes:]

        for comp in self.circuit.get_components():
            n1, n2, *extra_nodes = comp.connected_nodes

            # Condensatori: memorizzano la tensione ai loro capi
            if isinstance(comp, Capacitor):
                voltage_across_c = node_voltages[n1] - node_voltages[n2]
                comp.update_state(voltage_across_c, None) # La corrente viene ricalcolata nel build_mna_system
                                                          # o potresti calcolarla qui se vuoi memorizzarla.
            # Induttori: memorizzano la corrente attraverso di essi
            elif isinstance(comp, Inductor):
                inductor_aux_idx = self.inductor_aux_map[comp.component_id] - self.num_nodes
                current_through_l = aux_currents[inductor_aux_idx]
                comp.update_state(None, current_through_l) # La tensione viene ricalcolata in build_mna_system
                                                           # o potresti calcolarla qui.
            # LDR: la sua resistenza dinamica è già aggiornata dal livello di luce
            # (presupponendo che il livello di luce venga applicato esternamente prima della risoluzione del campione)
            # Non ha uno stato interno che evolve con la risoluzione del circuito, ma con un input esterno.
            elif isinstance(comp, LDR):
                pass # L'LDR aggiorna la sua resistenza internamente quando gli viene dato un livello di luce.
                     # Questo non è dipendente dalla risoluzione MNA, ma da un input esterno.


    def process_audio(self, input_signal_samples, input_node_name="Input", output_node_name="Output"):
        """
        Elabora un array di campioni audio attraverso il circuito simulato.

        Args:
            input_signal_samples (np.ndarray): Array di campioni del segnale di ingresso.
            input_node_name (str): Il nome del nodo al quale applicare il segnale di ingresso.
                                   Questo nodo sarà trattato come una sorgente di tensione ideale rispetto a GND.
            output_node_name (str): Il nome del nodo dal quale leggere il segnale di uscita.

        Returns:
            np.ndarray: Array dei campioni del segnale di uscita.
        """
        output_signal = np.zeros_like(input_signal_samples, dtype=float)

        # Inizializza lo stato del solutore per l'inizio della simulazione
        self.x_prev = np.zeros(self.total_vars)
        self.x_prev[self.circuit.get_ground_node_id()] = 0.0 # GND è sempre 0V

        # Ottieni gli ID dei nodi di input e output
        input_node_id = self.circuit.get_node_id(input_node_name)
        output_node_id = self.circuit.get_node_id(output_node_name)

        if input_node_id == self.circuit.get_ground_node_id():
            raise ValueError(f"Il nodo di input '{input_node_name}' non può essere GND.")
        if output_node_id == self.circuit.get_ground_node_id():
            print(f"Attenzione: Il nodo di output '{output_node_name}' è GND. L'uscita sarà 0.")

        print(f"\nInizio processamento audio per {len(input_signal_samples)} campioni...")
        for i, sample_value in enumerate(input_signal_samples):
            # Per applicare l'input, modifichiamo il guess iniziale forzando il nodo di input
            # Questo è un trucco comune per le sorgenti di tensione ideali.
            # Per un solutore più robusto, si aggiungerebbe una variabile ausiliaria per la corrente della sorgente.
            self.x_prev[input_node_id] = sample_value # Forza il valore del nodo di input

            # Risolvi il circuito per il campione corrente
            resolved_x = self.solve_sample(sample_value) # solve_sample userà l'input_signal_value per l'MNA

            # Estrai la tensione del nodo di output
            output_signal[i] = resolved_x[output_node_id]

            if (i + 1) % (self.sample_rate // 10) == 0: # Stampa avanzamento ogni 0.1s
                print(f"Elaborati {i+1}/{len(input_signal_samples)} campioni...", end='\r')
        print("\nProcessamento audio completato.")

        return output_signal


# Esempio di utilizzo (per testare la creazione delle istanze e i print)
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from components.resistor import Resistor
    from components.capacitor import Capacitor
    from components.diode import Diode
    from components.inductor import Inductor # Per testare le aux vars

    print("\n--- Test della classe CircuitSolver ---")
    my_test_circuit = Circuit("Filtro RC con Diodo")

    # Aggiungi nodi
    my_test_circuit.add_node("Input")
    my_test_circuit.add_node("Node1")
    my_test_circuit.add_node("Output")
    # GND è già aggiunto dal costruttore del Circuito

    # Aggiungi componenti
    # Resistor tra Input e Node1
    my_test_circuit.add_component(Resistor(1000), "Input", "Node1")
    # Diode tra Node1 e Output
    my_test_circuit.add_component(Diode(saturation_current=1e-12), "Node1", "Output")
    # Capacitor tra Output e GND
    my_test_circuit.add_component(Capacitor(1e-7, sample_rate=DEFAULT_SAMPLE_RATE), "Output", "GND")
    # Aggiungiamo un induttore per testare le variabili ausiliarie
    my_test_circuit.add_component(Inductor(10e-3, sample_rate=DEFAULT_SAMPLE_RATE), "Node1", "GND")

    print(my_test_circuit)

    # Inizializza il solutore con il circuito
    solver = CircuitSolver(my_test_circuit, sample_rate=DEFAULT_SAMPLE_RATE)

    # Genera un segnale di ingresso di test (onda sinusoidale)
    duration_s = 0.1 # secondi
    num_samples = int(DEFAULT_SAMPLE_RATE * duration_s)
    time_vector = np.linspace(0, duration_s, num_samples, endpoint=False)
    # Segnale di ingresso con ampiezza che vada oltre la caduta del diodo
    input_test_signal = 1.0 * np.sin(2 * np.pi * 100 * time_vector) # 1Vpk a 100Hz

    # Processa il segnale audio (il solutore MNA non è ancora realmente funzionante qui)
    output_test_signal = solver.process_audio(input_test_signal, input_node_name="Input", output_node_name="Output")

    # Plot dei risultati
    plt.figure(figsize=(12, 6))
    plt.plot(time_vector, input_test_signal, label='Segnale di Ingresso (V)')
    plt.plot(time_vector, output_test_signal, label='Segnale di Uscita (V) [Placeholder]')
    plt.title('Simulazione del Circuito (MNA - Logica MNA Implementata, Modelli Dettagliati da completare)')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Tensione (V)')
    plt.legend()
    plt.grid(True)
    plt.show()

    print("\nTest della classe CircuitSolver completato.")
    print("Nota: L'output mostra ancora un comportamento piatto o non realistico perché i metodi")
    print("di conduttanza/transconduttanza nei singoli componenti DEVONO essere implementati per")
    print("fornire le derivate analitiche necessarie alla Jacobiana di Newton-Raphson.")
    print("In assenza, viene usato un fallback numerico, ma l'accuratezza è compromessa.")
