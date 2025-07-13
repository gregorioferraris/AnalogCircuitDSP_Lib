import numpy as np
from utils.helpers import numerical_jacobian # Assicurati che questo percorso sia corretto e la funzione esista

class JFET:
    """
    Modello semplificato di JFET (Junction Field-Effect Transistor) N-channel.
    """
    def __init__(self, name: str, nodes: dict,
                 Idss: float = 0.01, Vp: float = -2.0, lambda_val: float = 0.05):
        """
        Args:
            name (str): Nome del componente.
            nodes (dict): Dizionario dei nodi: {'drain': 'node_D', 'gate': 'node_G', 'source': 'node_S'}.
            Idss (float): Corrente di Drain-Source a Vgs=0 e Vds > |Vp| (Ampere).
            Vp (float): Tensione di Pinch-Off (Volt). Per N-channel è negativo.
            lambda_val (float): Parametro di modulazione di canale (per Id vs Vds in saturazione).
        """
        if Idss <= 0:
            raise ValueError("Idss deve essere positivo.")
        if Vp >= 0: # Per JFET N-channel, Vp deve essere negativo
            raise ValueError("Vp (tensione di pinch-off) deve essere negativo per JFET N-channel.")
        if lambda_val < 0:
            raise ValueError("Lambda (modulazione di canale) deve essere non negativo.")

        self.name = name
        self.nodes = nodes
        self.Idss = float(Idss)
        self.Vp = float(Vp)
        self.lambda_val = float(lambda_val)

        required_nodes = ['drain', 'gate', 'source']
        if not all(node in self.nodes for node in required_nodes):
            raise ValueError(f"JFET '{name}': Nodi richiesti (drain, gate, source) non specificati correttamente.")

        # Memorizzazione dello stato precedente delle tensioni per debug o modelli complessi
        self.last_Vgs = 0.0
        self.last_Vds = 0.0

    def calculate_drain_current(self, v_gs: float, v_ds: float) -> float:
        """
        Calcola la corrente di Drain (Id) del JFET.
        Utilizza un modello semplificato che considera regione di triodo/saturazione.
        """
        self.last_Vgs = v_gs
        self.last_Vds = v_ds

        # Gestione della giunzione gate-source polarizzata in avanti
        if v_gs >= 0.0: # Per un JFET a canale N, Vgs positiva polarizza in avanti la giunzione
            # Questa è una semplificazione. Un modello completo considererebbe la corrente di gate (Ig)
            # e l'impatto sulla Id. Per ora, limitiamo la Id se Vgs è positiva.
            # Potresti voler implementare un modello di diodo qui per Ig.
            # Per il momento, assumiamo che non ci sia conduzione significativa del canale se Vgs è troppo alta e positiva.
            # In genere un JFET non è usato con Vgs > 0, se non per valori molto piccoli.
            # Una Id = 0.0 è troppo drastica per Vgs > 0.
            # Un approccio migliore potrebbe essere limitare Vgs per il calcolo della Id a un piccolo valore positivo.
            # Tuttavia, per ora, per semplicità e per evitare comportamenti complessi non modellati, manteniamo 0.
            return 0.0

        if v_gs <= self.Vp:
            # Regione di Cut-off (Vgs < Vp)
            return 0.0
        elif v_ds >= (v_gs - self.Vp):
            # Regione di Saturazione (Vds >= Vgs - Vp)
            # Id = Idss * (1 - Vgs/Vp)^2 * (1 + lambda * Vds)
            term_vgs = (1.0 - v_gs / self.Vp)
            if term_vgs < 0: term_vgs = 0 # Clamping per robustezza
            id_saturation = self.Idss * (term_vgs ** 2) * (1.0 + self.lambda_val * v_ds)
            return max(0.0, id_saturation) # La corrente non può essere negativa
        else:
            # Regione di Triodo (Ohmica) (Vds < Vgs - Vp)
            # Id = Idss * [2(1 - Vgs/Vp)Vds - Vds^2]/Vp^2 * (1 + lambda * Vds) - Lambda è spesso ignorato in triodo
            # Ho aggiunto (1 + lambda * Vds) per coerenza, ma molti modelli lo omettono in triodo.
            V_effective_pinchoff = v_gs - self.Vp
            # Evitare divisione per zero o problemi con Vp vicino a zero se ci fossero casi limite
            if self.Vp == 0: return 0.0 # O gestire errore
            
            # Assicurati che V_effective_pinchoff sia positivo per la formula.
            # Per Vds < Vgs - Vp, V_effective_pinchoff > Vds.
            # La formula del triodo:
            # Id = Beta * [2 * (Vgs - Vth) * Vds - Vds^2]
            # dove Beta = Idss / Vp^2
            Beta = self.Idss / (self.Vp**2)
            id_triode = Beta * (2 * (v_gs - self.Vp) * v_ds - v_ds**2) * (1.0 + self.lambda_val * v_ds)
            return max(0.0, id_triode)

    def calculate_transconductance(self, v_gs: float, v_ds: float) -> float:
        """
        Calcola la transconduttanza (gm = d(Id)/d(Vgs)) del JFET usando derivazione numerica.
        """
        # Avvolgi la funzione in una lambda che accetta solo v_gs
        func_to_derive = lambda current_v_gs: self.calculate_drain_current(current_v_gs, v_ds)
        
        # Aggiungi un piccolo offset per la derivata numerica
        epsilon = 1e-6 
        return numerical_jacobian(func_to_derive, v_gs, h=epsilon)

    def calculate_output_conductance(self, v_gs: float, v_ds: float) -> float:
        """
        Calcola la conduttanza di output (gds = d(Id)/d(Vds)) del JFET usando derivazione numerica.
        """
        # Avvolgi la funzione in una lambda che accetta solo v_ds
        func_to_derive = lambda current_v_ds: self.calculate_drain_current(v_gs, current_v_ds)

        # Aggiungi un piccolo offset per la derivata numerica
        epsilon = 1e-6
        return numerical_jacobian(func_to_derive, v_ds, h=epsilon)

    def calculate_jacobian_elements(self, V_drain: float, V_gate: float, V_source: float) -> np.ndarray:
        """
        Calcola gli elementi della matrice Jacobiana per il JFET.
        La matrice è 3x3, con righe/colonne corrispondenti a [Drain, Gate, Source].
        
        Parametri:
            V_drain (float): Tensione assoluta al Drain.
            V_gate (float): Tensione assoluta al Gate.
            V_source (float): Tensione assoluta al Source.
            
        Returns:
            np.ndarray: Matrice Jacobiana 3x3.
        """
        # Calcola le tensioni differenziali per i metodi interni del JFET
        v_gs_current = V_gate - V_source
        v_ds_current = V_drain - V_source

        # Se in cutoff o Vgs positiva, le derivate sono zero.
        if v_gs_current < self.Vp or v_gs_current >= 0.0:
            return np.zeros((3, 3))

        # Calcola le conduttanze numericamente
        gm = self.calculate_transconductance(v_gs_current, v_ds_current)
        gds = self.calculate_output_conductance(v_gs_current, v_ds_current)

        # Matrice Jacobiana 3x3: [Drain, Gate, Source]
        jacobian_matrix = np.zeros((3, 3))

        # Corrente di Drain (Id) entra nel Drain. F_drain += Id.
        # Derivate di Id rispetto ai nodi:
        # d(Id)/d(V_drain) = gds
        # d(Id)/d(V_gate) = gm
        # d(Id)/d(V_source) = -gm - gds (poiché Vgs e Vds dipendono da Vsource)
        jacobian_matrix[0, 0] = gds      # d(F_drain)/d(V_drain)
        jacobian_matrix[0, 1] = gm       # d(F_drain)/d(V_gate)
        jacobian_matrix[0, 2] = -gm - gds # d(F_drain)/d(V_source)

        # Corrente di Gate (Ig) è assunta zero per un JFET ideale. F_gate += Ig (0).
        # Quindi la riga della Jacobiana per il Gate è tutta zero.
        # jacobian_matrix[1, :] = 0.0

        # Corrente di Source (Is) esce dal Source. F_source -= Is.
        # Poiché Id + Ig + Is = 0 e Ig = 0, allora Is = -Id.
        # Quindi, d(F_source)/dVx = -d(Is)/dVx = -(-d(Id)/dVx) = d(Id)/dVx.
        # Le derivate del nodo Source sono le stesse del nodo Drain, ma con segno invertito se consideriamo Is = -(Id + Ig).
        # No, sono uguali a quelle del Drain perché F_source = -Is e Is = -Id
        # Quindi d(F_source)/dVx = d(Id)/dVx.
        jacobian_matrix[2, 0] = gds      # d(F_source)/d(V_drain)
        jacobian_matrix[2, 1] = gm       # d(F_source)/d(V_gate)
        jacobian_matrix[2, 2] = -gm - gds # d(F_source)/d(V_source)

        return jacobian_matrix

    def set_nodes(self, drain_node: str, gate_node: str, source_node: str):
        """
        Imposta i nomi dei nodi a cui il JFET è collegato nel circuito.
        """
        self.nodes['drain'] = drain_node
        self.nodes['gate'] = gate_node
        self.nodes['source'] = source_node

    def get_node_names(self) -> dict:
        """
        Restituisce i nomi dei nodi a cui il JFET è collegato.
        """
        return self.nodes

    def __str__(self):
        return f"JFET(name='{self.name}', Idss={self.Idss:.1e}A, Vp={self.Vp:.1f}V, lambda={self.lambda_val:.2f})"

# Esempio di utilizzo (per testare la classe JFET in isolamento)
if __name__ == "__main__":
    # Mini-implementazione di numerical_jacobian per il test locale
    def numerical_jacobian(func, x0, h=1e-6):
        return (func(x0 + h) - func(x0 - h)) / (2 * h)

    # Parametri tipici per un JFET a canale N (es. J201)
    jfet_params = {
        "Idss": 0.003,  # 3 mA
        "Vp": -0.8,     # -0.8 V
        "lambda_val": 0.01 # Piccolo effetto di modulazione di canale
    }

    my_jfet = JFET("My_J201_JFET", 
                   nodes={'drain': 'D1', 'gate': 'G1', 'source': 'S1'},
                   **jfet_params)

    print(f"JFET: {my_jfet.name}")
    print(f"Parameters: Idss={my_jfet.Idss*1000:.3f}mA, Vp={my_jfet.Vp}V, lambda={my_jfet.lambda_val}")

    # Esempi di calcolo della corrente di drain
    V_gs_test = -0.4
    V_ds_test = 5.0
    Id = my_jfet.calculate_drain_current(V_gs_test, V_ds_test)
    print(f"\nVgs = {V_gs_test}V, Vds = {V_ds_test}V -> Id = {Id*1000:.3f} mA (Saturazione)")

    V_gs_test = -1.0 # In cutoff
    V_ds_test = 5.0
    Id = my_jfet.calculate_drain_current(V_gs_test, V_ds_test)
    print(f"Vgs = {V_gs_test}V (cutoff) -> Id = {Id*1000:.3f} mA")
    
    V_gs_test = 0.0 # Vgs = 0, Id dovrebbe essere Idss
    V_ds_test = 5.0
    Id = my_jfet.calculate_drain_current(V_gs_test, V_ds_test)
    print(f"Vgs = {V_gs_test}V (Vgs=0) -> Id = {Id*1000:.3f} mA (dovrebbe essere Idss * (1 + lambda * Vds))")

    V_gs_test = -0.4 # Triodo region
    V_ds_test = 0.1 # Vds < Vgs - Vp (0.1 < -0.4 - (-0.8) = 0.4)
    Id = my_jfet.calculate_drain_current(V_gs_test, V_ds_test)
    print(f"Vgs = {V_gs_test}V, Vds = {V_ds_test}V (Triodo) -> Id = {Id*1000:.3f} mA")

    # Esempio di calcolo della Jacobiana
    V_drain_curr = 5.0
    V_gate_curr = -0.4
    V_source_curr = 0.0
    
    jacobian = my_jfet.calculate_jacobian_elements(V_drain_curr, V_gate_curr, V_source_curr)
    print(f"\nJacobian Matrix for Vd={V_drain_curr}V, Vg={V_gate_curr}V, Vs={V_source_curr}V:")
    print(jacobian * 1000) # Per visualizzare in mS
    print("(Valori in mS, Drain, Gate, Source)")

    V_drain_curr = 5.0
    V_gate_curr = -1.0 # In cutoff
    V_source_curr = 0.0
    jacobian = my_jfet.calculate_jacobian_elements(V_drain_curr, V_gate_curr, V_source_curr)
    print(f"\nJacobian Matrix (cutoff) for Vd={V_drain_curr}V, Vg={V_gate_curr}V, Vs={V_source_curr}V:")
    print(jacobian * 1000)
