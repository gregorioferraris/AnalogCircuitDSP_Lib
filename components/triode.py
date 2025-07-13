import numpy as np

class Triode:
    """
    Modello di triodo per simulazioni di circuiti elettronici audio.
    Basato su un modello comune che deriva dalla caratteristica di transconduttanza.
    """
    def __init__(self, name: str, nodes: dict,
                 mu: float = 100.0, # Fattore di amplificazione (gain factor)
                 K: float = 0.001,  # Costante di emissione (emission constant)
                 x: float = 1.0,    # Esponente, spesso 1.0 (se lineare) o 1.5 per modelli ideali
                 Vg_cutoff: float = -1.5): # Tensione di cutoff della griglia (grid cutoff voltage)
        
        self.name = name
        # I nodi dovrebbero essere un dizionario: {'anode': 'node_A', 'grid': 'node_G', 'cathode': 'node_K'}
        self.nodes = nodes
        
        self.mu = mu
        self.K = K
        self.x = x
        self.Vg_cutoff = Vg_cutoff

        # Assicurati che i nodi necessari siano presenti
        required_nodes = ['anode', 'grid', 'cathode']
        if not all(node in self.nodes for node in required_nodes):
            raise ValueError(f"Triode '{name}': Nodi richiesti (anode, grid, cathode) non specificati correttamente.")

        # Inizializza lo stato (non richiesto direttamente dal solutore, ma utile per modelli complessi)
        self.last_V_gk = 0.0
        self.last_V_ak = 0.0

    def calculate_plate_current(self, V_gk: float, V_ak: float) -> float:
        """
        Calcola la corrente di placca (anodo) del triodo.
        
        Args:
            V_gk (float): Tensione griglia-catodo.
            V_ak (float): Tensione anodo-catodo.
            
        Returns:
            float: Corrente di placca (in Ampere).
        """
        self.last_V_gk = V_gk # Memorizza per debug o futuro
        self.last_V_ak = V_ak # Memorizza per debug o futuro

        # Modello di Koren/Log per triodi (semplificato)
        # Ef = V_gk + V_ak / mu
        # Se Ef < Vg_cutoff, la corrente è 0 (cutoff)
        # Ia = K * (Ef - Vg_cutoff)^x
        
        # Una versione semplificata che cattura il comportamento non lineare:
        # Tensione effettiva alla griglia
        V_eff = V_gk + V_ak / self.mu

        # Applica la condizione di cutoff
        if V_eff < self.Vg_cutoff:
            return 0.0
        
        # Calcola la corrente di placca. Applica un piccolo offset per evitare log(0) o pow(negativo) se `x` non intero.
        # np.maximum(0, ...) assicura che il termine sotto la potenza non sia negativo.
        Ia = self.K * (np.maximum(0, V_eff - self.Vg_cutoff))**self.x
        
        return Ia

    def calculate_jacobian_elements(self, V_anode: float, V_grid: float, V_cathode: float) -> np.ndarray:
        """
        Calcola gli elementi della matrice Jacobiana per il triodo.
        Questa matrice rappresenta le derivate parziali della corrente di placca (Ia)
        rispetto alle tensioni ai suoi nodi (Anode, Grid, Cathode).
        
        Derivate necessarie:
        d(Ia)/d(V_anode)
        d(Ia)/d(V_grid)
        d(Ia)/d(V_cathode)
        
        La matrice Jacobiana per il triodo è una 3x3, con righe/colonne corrispondenti a [Anode, Grid, Cathode].
        Tuttavia, poiché la corrente di placca fluisce solo tra anodo e catodo
        (e la corrente di griglia è assunta ideale = 0),
        la Jacobiana si riempirà solo per le righe/colonne pertinenti.
        
        La matrice J_comp sarà:
        [[ dIa/dVa, dIa/dVg, dIa/dVk  ]  (Contributo nodo Anodo)
         [ 0,       0,       0      ]  (Contributo nodo Grid - se corrente griglia = 0)
         [ dIk/dVa, dIk/dVg, dIk/dVk  ]] (Contributo nodo Cathode)
         
        Dato che Ia = -Ik (corrente anodo è opposta alla corrente catodo), dIk/dVx = -dIa/dVx.
        
        V_gk = V_grid - V_cathode
        V_ak = V_anode - V_cathode
        V_eff = V_gk + V_ak / mu = V_grid - V_cathode + (V_anode - V_cathode) / mu
              = V_grid + V_anode / mu - V_cathode * (1 + 1/mu)

        Ia = K * (max(0, V_eff - Vg_cutoff))^x
        
        Deriviamo Ia rispetto a Va, Vg, Vk:
        
        d(Ia)/d(V_anode) = K * x * (V_eff - Vg_cutoff)^(x-1) * d(V_eff)/d(V_anode)
                         = K * x * (V_eff - Vg_cutoff)^(x-1) * (1/mu)
        Questa è la transconduttanza di placca (plate conductance), spesso 1/rp.
        
        d(Ia)/d(V_grid)  = K * x * (V_eff - Vg_cutoff)^(x-1) * d(V_eff)/d(V_grid)
                         = K * x * (V_eff - Vg_cutoff)^(x-1) * (1)
        Questa è la transconduttanza (transconductance) Gm.
        
        d(Ia)/d(V_cathode) = K * x * (V_eff - Vg_cutoff)^(x-1) * d(V_eff)/d(V_cathode)
                           = K * x * (V_eff - Vg_cutoff)^(x-1) * (-(1 + 1/mu))
        
        """
        # Ricalcola V_gk e V_ak basandosi sulle tensioni attuali passate
        V_gk_current = V_grid - V_cathode
        V_ak_current = V_anode - V_cathode
        
        V_eff = V_gk_current + V_ak_current / self.mu

        # Se il triodo è in cutoff, le derivate sono 0.
        if V_eff < self.Vg_cutoff:
            return np.zeros((3, 3)) # Matrice 3x3 di zeri per [Anode, Grid, Cathode]

        # Termine comune alle derivate
        common_term = self.K * self.x * (np.maximum(0, V_eff - self.Vg_cutoff))**(self.x - 1)
        
        # Derivate parziali di Ia
        dIa_dVa = common_term * (1.0 / self.mu) # Transconduttanza di placca (1/rp)
        dIa_dVg = common_term * 1.0              # Transconduttanza (Gm)
        dIa_dVc = common_term * (-(1.0 + 1.0 / self.mu)) # Derivata rispetto al catodo

        # La Jacobiana per un componente non lineare con N nodi coinvolti
        # è una matrice NxN dove J_ij = d(Corrente_che_esce_da_i) / d(V_j)
        # Oppure, per il lato destro dell'MNA (F(x)), è la derivata di F_i rispetto a V_j
        # F_i rappresenta la somma delle correnti che ENTRANO nel nodo i.
        # Se Ia entra nell'anodo e esce dal catodo, allora F_anode += Ia, F_cathode -= Ia
        # d(F_anode)/dVx = d(Ia)/dVx
        # d(F_cathode)/dVx = -d(Ia)/dVx

        jacobian_matrix = np.zeros((3, 3)) # Indici: [Anode, Grid, Cathode]

        # Contributi per la riga dell'Anodo (derivata della corrente del nodo Anodo)
        # d(F_anode)/d(V_anode) = d(Ia)/d(V_anode)
        jacobian_matrix[0, 0] = dIa_dVa 
        # d(F_anode)/d(V_grid) = d(Ia)/d(V_grid)
        jacobian_matrix[0, 1] = dIa_dVg
        # d(F_anode)/d(V_cathode) = d(Ia)/d(V_cathode)
        jacobian_matrix[0, 2] = dIa_dVc

        # Contributi per la riga del Catodo (derivata della corrente del nodo Catodo)
        # Poiché Ik = -Ia, d(F_cathode)/dVx = d(-Ia)/dVx = -d(Ia)/dVx
        jacobian_matrix[2, 0] = -dIa_dVa
        jacobian_matrix[2, 1] = -dIa_dVg
        jacobian_matrix[2, 2] = -dIa_dVc

        # La riga della Griglia è zero se non consideriamo corrente di griglia.
        # Se si volesse modellare la corrente di griglia (es. per bias in classe A2),
        # si dovrebbe aggiungere una funzione calculate_grid_current e le sue derivate.

        return jacobian_matrix

    def set_nodes(self, anode_node: str, grid_node: str, cathode_node: str):
        """
        Imposta i nomi dei nodi a cui il triodo è collegato nel circuito.
        """
        self.nodes['anode'] = anode_node
        self.nodes['grid'] = grid_node
        self.nodes['cathode'] = cathode_node

    def get_node_names(self) -> dict:
        """
        Restituisce i nomi dei nodi a cui il triodo è collegato.
        """
        return self.nodes

# Esempio di utilizzo (per testare la classe Triode in isolamento)
if __name__ == "__main__":
    # Parametri tipici per una 12AX7
    triode_params = {
        "mu": 100.0,
        "K": 80.0e-6, # Da calibrare per il modello esatto
        "x": 1.5,     # Valore comune per esponente
        "Vg_cutoff": -1.2
    }

    # Creazione di un'istanza di triodo
    my_triode = Triode("My_12AX7_Triode", 
                       nodes={'anode': 'ANODE_NODE', 'grid': 'GRID_NODE', 'cathode': 'CATHODE_NODE'},
                       **triode_params)

    print(f"Triode: {my_triode.name}")
    print(f"Parameters: mu={my_triode.mu}, K={my_triode.K}, x={my_triode.x}, Vg_cutoff={my_triode.Vg_cutoff}")

    # Esempi di calcolo della corrente di placca
    V_gk_test = -1.0
    V_ak_test = 150.0
    Ia = my_triode.calculate_plate_current(V_gk_test, V_ak_test)
    print(f"\nV_gk = {V_gk_test}V, V_ak = {V_ak_test}V -> Ia = {Ia*1000:.3f} mA")

    V_gk_test = -3.0 # In cutoff
    V_ak_test = 150.0
    Ia = my_triode.calculate_plate_current(V_gk_test, V_ak_test)
    print(f"V_gk = {V_gk_test}V, V_ak = {V_ak_test}V (cutoff) -> Ia = {Ia*1000:.3f} mA")

    V_gk_test = 0.0
    V_ak_test = 200.0
    Ia = my_triode.calculate_plate_current(V_gk_test, V_ak_test)
    print(f"V_gk = {V_gk_test}V, V_ak = {V_ak_test}V -> Ia = {Ia*1000:.3f} mA")
    
    # Esempio di calcolo della Jacobiana
    V_anode_curr = 150.0
    V_grid_curr = -1.0
    V_cathode_curr = 0.0
    
    jacobian = my_triode.calculate_jacobian_elements(V_anode_curr, V_grid_curr, V_cathode_curr)
    print(f"\nJacobian Matrix for V_anode={V_anode_curr}V, V_grid={V_grid_curr}V, V_cathode={V_cathode_curr}V:")
    print(jacobian * 1000) # Per visualizzare in mS
    print("(Valori in mS, Anodo, Griglia, Catodo)")

    V_anode_curr = 150.0
    V_grid_curr = -3.0 # In cutoff
    V_cathode_curr = 0.0
    jacobian = my_triode.calculate_jacobian_elements(V_anode_curr, V_grid_curr, V_cathode_curr)
    print(f"\nJacobian Matrix (cutoff) for V_anode={V_anode_curr}V, V_grid={V_grid_curr}V, V_cathode={V_cathode_curr}V:")
    print(jacobian * 1000)
