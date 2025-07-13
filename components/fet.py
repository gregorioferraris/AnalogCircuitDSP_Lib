# components/fet.py

import numpy as np

class FET:
    """
    Modello semplificato di un JFET (Junction Field-Effect Transistor)
    per l'uso come Resistore Controllato in Tensione (VCR) nella regione ohmica.
    Il FET agisce come una resistenza controllata dalla tensione Vgs.

    Parametri tipici per un JFET come 2N5457, 2N5484, J201.
    """
    def __init__(self, name, vt=-2.0, kp=0.005, channel_resistance_at_zero_vgs=100.0):
        """
        Inizializza il FET.
        Args:
            name (str): Nome del componente.
            vt (float): Tensione di soglia (pinch-off voltage) Vp.
                        La tensione Vgs al di sotto della quale il FET si spegne (n-channel) o si accende (p-channel).
                        Per JFET n-channel è negativa.
            kp (float): Coefficiente di transconduttanza (conducibilità) in A/V^2.
                        Determina la pendenza della curva Id vs Vgs.
            channel_resistance_at_zero_vgs (float): Resistenza del canale D-S quando Vgs = 0V (regione ohmica).
                                                   Usata per calibrare il comportamento.
        """
        self.name = name
        self.vt = float(vt)  # Vp, tensione di pinch-off (gate-source)
        self.kp = float(kp)  # Transconduttanza (beta)
        self.Rds_on_at_zero_vgs = float(channel_resistance_at_zero_vgs) # R_on quando Vgs = 0

        # Parametri interni per il calcolo della Jacobiana
        self._vgs_current = 0.0
        self._vds_current = 0.0

        # Mappatura dei nodi (Source, Gate, Drain)
        self.nodes = {} # Mappa il nome del pin al Node ID (gestito dal circuito)

    def set_nodes(self, source_node_id, gate_node_id, drain_node_id):
        """Assegna gli ID dei nodi ai pin del FET."""
        self.nodes['source'] = source_node_id
        self.nodes['gate'] = gate_node_id
        self.nodes['drain'] = drain_node_id

    def calculate_currents(self, v_source, v_gate, v_drain):
        """
        Calcola le correnti di Drain-Source (Ids) e Gate (Ig)
        in base alle tensioni ai nodi.
        Ci concentriamo sulla regione ohmica per la compressione.

        Per un JFET n-channel:
        Ids = kp * ((Vgs - Vt) * Vds - Vds^2 / 2)  (regione ohmica, Vds < Vgs - Vt)
        Assumiamo Vgs_eff = Vgs - Vt
        Ids = kp * (Vgs_eff * Vds - Vds^2 / 2)

        Per l'MNA, è più semplice modellare la resistenza Drain-Source (Rds) e
        usare Ids = Vds / Rds.
        Rds = R_on / (1 - Vgs/Vp)^2 (modello approssimato per JFET)
        Dove R_on è la resistenza con Vgs=0.
        Questo modello è più adatto per il VCR.
        """
        vgs = v_gate - v_source
        vds = v_drain - v_source

        # Assicuriamoci che Vgs non superi 0V per JFET N-channel (tipicamente Vgs <= 0)
        # E che Vgs sia sopra la tensione di pinch-off (Vt)
        vgs_clamped = np.clip(vgs, self.vt, 0.0) # Vgs deve essere tra Vp e 0V

        # Calcolo della resistenza Drain-Source (Rds)
        # Modello semplificato di resistenza variabile per regione ohmica:
        if vgs_clamped <= self.vt: # Il FET è spento (alta impedenza)
            rds = 1e12 # Molto alta resistenza
        else:
            # Formula più robusta per Rds_on in base a Vgs:
            # Rds = Rds_on_at_zero_vgs / (1 - Vgs/Vp)^2
            # Per JFET, Vp è la tensione di pinch-off (Vt), spesso negativa.
            # Assumiamo Vgs/Vp è un rapporto, dove Vp < 0.
            # (1 - Vgs/Vt)
            ratio = (1.0 - vgs_clamped / self.vt) # Questo rapporto va da 0 a 1 (quando Vgs=0)
            if ratio < 0.01: # Evita divisioni per zero o valori molto piccoli
                ratio = 0.01
            rds = self.Rds_on_at_zero_vgs / (ratio * ratio)
            
            # Limita la resistenza massima per evitare problemi numerici
            rds = np.clip(rds, self.Rds_on_at_zero_vgs, 1e9) # Max 1 GigaOhm

        # Corrente Drain-Source
        ids = vds / rds

        # Corrente di Gate (tipicamente trascurabile per JFET ideali, ma può avere una piccola perdita)
        # Per MNA, la modellazione precisa include una resistenza Rin per il gate.
        ig = 0.0 # Per un modello ideale, la corrente di gate è zero.

        self._vgs_current = vgs # Memorizza per la Jacobiana
        self._vds_current = ids # Memorizza per la Jacobiana

        return ids, ig, rds # Restituisce anche Rds per riferimenti esterni

    def calculate_jacobian(self, v_source, v_gate, v_drain):
        """
        Calcola gli elementi della matrice Jacobiana per il FET.
        Si tratta delle derivate parziali delle correnti rispetto alle tensioni ai nodi.
        J_ij = d(I_i) / d(V_j)
        Dove I_i è la corrente che fluisce dal nodo i, e V_j è la tensione al nodo j.

        Per un FET come VCR, stiamo modellando principalmente la conduttanza Gds = Ids / Vds.
        Gds = 1 / Rds.
        La derivata di Ids rispetto a Vds è 1/Rds (se Rds è costante rispetto a Vds)
        La derivata di Ids rispetto a Vgs è d(Ids)/d(Vgs).

        Questo è il punto più complesso per la non linearità.
        Utilizzeremo una derivata numerica o un modello semplificato.
        """
        vgs = v_gate - v_source
        vds = v_drain - v_source

        # Conduttanza Drain-Source (Gds)
        # Ids = Vds * Gds(Vgs)
        # d(Ids)/d(Vds) = Gds(Vgs)
        # d(Ids)/d(Vgs) = Vds * d(Gds)/d(Vgs)

        # Usiamo il metodo numerico per la Jacobiana per questo componente non lineare.
        # Questa funzione dovrebbe restituire una sottomatrice Jacobiana 3x3
        # relativa a (Drain, Gate, Source) o ai nodi interni del circuito.
        
        # J = [ [ dId/dVd, dId/dVg, dId/dVs ],
        #       [ dIg/dVd, dIg/dVg, dIg/dVs ],
        #       [ dIs/dVd, dIs/dVg, dIs/dVs ] ]
        # Con Is = -(Id + Ig)

        # Data la complessità, è comune approssimare d(Ids)/d(Vgs) o usarne un modello semplificato.
        # Per ora, restituiamo un modello base che usa la conduttanza inversa Rds.
        
        # Riutilizziamo la logica di calcolo delle correnti per ottenere rds_current
        ids, ig, rds_current = self.calculate_currents(v_source, v_gate, v_drain)
        
        # Matrice Jacobiana 3x3 per (Drain, Gate, Source)
        jacobian = np.zeros((3, 3))

        # Nodo Drain (0)
        # Corrente I_drain = Ids (che entra nel drain se Vds > 0)
        # d(I_drain)/d(V_drain) = 1 / Rds
        jacobian[0, 0] = 1.0 / rds_current # d(Ids)/d(Vd)
        jacobian[0, 1] = 0.0                # d(Ids)/dVg (approssimato, più complesso)
        jacobian[0, 2] = -1.0 / rds_current # d(Ids)/dVs

        # Nodo Gate (1)
        # Corrente I_gate = 0 (ideale) -> tutte le derivate sono 0
        # Se modellassimo una resistenza di ingresso del gate, sarebbe diversa.
        jacobian[1, 0] = 0.0
        jacobian[1, 1] = 0.0
        jacobian[1, 2] = 0.0

        # Nodo Source (2)
        # Corrente I_source = -Ids (che esce dal source)
        jacobian[2, 0] = -1.0 / rds_current # d(-Ids)/dVd
        jacobian[2, 1] = 0.0                 # d(-Ids)/dVg (approssimato)
        jacobian[2, 2] = 1.0 / rds_current  # d(-Ids)/dVs

        return jacobian

    def __str__(self):
        return f"FET({self.name}, Vp={self.vt:.2f}V, Kp={self.kp:.3e}, Rds_on_at_0Vgs={self.Rds_on_at_zero_vgs:.1f} Ohm)"
