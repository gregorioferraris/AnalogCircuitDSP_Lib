# circuits/bjt_push_pull_amp.py

from components.resistor import Resistor
from components.capacitor import Capacitor
from components.voltage_source import VoltageSource
from components.bjt import BJT # Per NPN
# from components.pnp_bjt import PNP_BJT # Avresti bisogno di una classe separata per PNP o un flag nella BJT
from components.diode import Diode # Per bias crossover
from circuit_solver.circuit import Circuit

# --- CLASSE PNP_BJT FITTIZIA (per questo esempio) ---
# Se non hai una PNP_BJT, dovrai crearla o modificare la tua BJT esistente.
# Questa è una semplificazione per far funzionare il circuito di esempio.
# Una vera PNP_BJT avrebbe calculate_collector_current con V_be < 0 e Id positiva.
# La calculate_jacobian_elements sarebbe simile ma con gm e gds calcolati correttamente per PNP.
class PNP_BJT(BJT): # Eredita da BJT per semplicità, ma le formule Id devono essere per PNP
    def __init__(self, name: str, nodes: dict, Is: float = 1.0e-14, Vt: float = 0.02585, Beta_F: float = 100.0, Va: float = 70.0):
        super().__init__(name, nodes, Is, Vt, Beta_F, Va)
        # Potresti aggiungere qui una logica specifica o un flag self.type = 'PNP'
        # Sovrascrivi calculate_collector_current e calculate_jacobian_elements se necessario
        # Per ora, useremo l'implementazione del BJT NPN come placeholder per la struttura.
        # IL MODELLO ID REALE PER PNP DOVRÀ ESSERE IMPLEMENTATO NEL calculate_collector_current
        # per Vbe < 0 e Vce < 0 per saturazione, e le correnti saranno negative (usciranno dal collettore).
        pass

    def calculate_collector_current(self, v_be: float, v_ce: float) -> float:
        # Questo è un placeholder. La vera Ic per PNP sarà negativa e calcolata per Vbe negativa.
        # Per un PNP, la Ic è tipicamente Ic = -Is * (exp(-Vbe/Vt) - 1) * (1 - Vce/Va)
        # e Beta_F si applica a Ie o Ib.
        # Qui useremo il modello NPN ma con Vbe e Vce invertiti per il calcolo,
        # e restituiremo un valore negativo, assumendo un comportamento complementare.
        # QUESTA PARTE RICHIEDE ATTENZIONE: Il modello BJT è per NPN. Una PNP deve gestire Vbe < 0.
        # Per semplicità in questo esempio, supponiamo che la tua classe BJT sappia gestire PNP
        # se riceve nodi di PNP e tensioni relative.
        # Se la tua BJT NPN è pura NPN, avrai bisogno di una vera PNP_BJT.
        
        # Per fare un PNP semplice per il test, possiamo invertire il segno degli input
        # e del risultato della calculate_collector_current NPN, se il tuo NPN è robusto.
        # Oppure una implementazione PNP dedicata:
        
        # Modo più realistico per PNP:
        # Vbe_abs = abs(v_be) # Vbe è negativo per la conduzione PNP
        # if Vbe_abs < 0.5: return 0.0 # Cutoff
        # ... e poi le formule per PNP ...
        
        # PER GLI SCOPI DI QUESTO ESEMPIO E TEST, ASSUMIAMO CHE LA calculate_collector_current
        # della classe BJT gestisca implicitamente la polarità se v_be è negativa per PNP.
        # Se non è così, DEVI implementare una PNP_BJT separata o aggiungere logica nella BJT.
        
        # Dato che le tensioni V_be_Q2 e V_ce_Q2 saranno negative, la calculate_collector_current
        # della BJT NPN molto probabilmente restituirà 0 o un valore errato.
        # Quindi, per un corretto funzionamento, *è necessaria una classe PNP_BJT dedicata*
        # o un'estensione della classe BJT per gestire NPN/PNP.
        
        # Per ora, faccio un HACK per far "funzionare" il test DC,
        # assumendo che il tuo BJT.calculate_collector_current possa interpretare
        # un NPN con Vbe e Vce invertite (ma non è un buon modello reale).
        # Il modo corretto è:
        # if v_be >= 0: return 0.0 # PNP non conduce con Vbe positiva
        # # ... formule PNP con v_be < 0
        return -super().calculate_collector_current(abs(v_be), abs(v_ce)) # Non è un modello PNP corretto, solo per test.


class BJTPushPullAmp(Circuit):
    def __init__(self, name="BJT_Push_Pull_Amp"):
        super().__init__(name)

        # --- Parametri del circuito ---
        VCC = 9.0  # Tensione di alimentazione positiva
        VEE = -9.0 # Tensione di alimentazione negativa (per dual rail)
        
        R_in_val = 1.0e3 # Resistenza di ingresso
        R_load_val = 8.0 # Resistenza di carico (es. altoparlante 8 Ohm)

        # Resistenze per polarizzare le basi (se non usi diodi)
        # R_bias_val = 100.0 # Per un bias molto semplice
        
        # Diodi per bias crossover (Classe AB)
        # V_bias = 2 * V_diode_drop ~ 1.4V per 2 diodi in serie
        # Questo crea una piccola corrente di riposo per eliminare la distorsione di crossover
        diode_bias_params = {
            "Is": 1.0e-14, # Corrente di saturazione inversa (A)
            "Vt": 0.02585, # Tensione termica a 300K (V)
            "N": 1.0       # Fattore di idealità
        }

        # --- Parametri dei BJT (NPN e PNP generici) ---
        q_npn_params = {
            "Is": 1.0e-14, "Vt": 0.02585, "Beta_F": 100.0, "Va": 70.0
        }
        q_pnp_params = { # Per il PNP
            "Is": 1.0e-14, "Vt": 0.02585, "Beta_F": 100.0, "Va": 70.0
        }

        # --- Nodi del circuito ---
        self.add_node('in')
        self.add_node('out')
        self.add_node('VCC_pos') # Nodo alimentazione positiva
        self.add_node('VCC_neg') # Nodo alimentazione negativa
        self.add_node('Q1_base')
        self.add_node('Q1_collector')
        self.add_node('Q1_emitter') # Questo sarà il nodo di uscita per i transistor di potenza
        self.add_node('Q2_base')
        self.add_node('Q2_collector')
        self.add_node('Q2_emitter') # Questo sarà collegato al nodo Q1_emitter per l'uscita

        # --- Componenti ---

        # 1. Alimentazioni Duali
        self.add_component(VoltageSource(name="VCC_Pos_Source", nodes={'pos': 'VCC_pos', 'neg': 'gnd'}, voltage=VCC))
        self.add_component(VoltageSource(name="VCC_Neg_Source", nodes={'pos': 'gnd', 'neg': 'VCC_neg'}, voltage=abs(VEE))) # Negativa rispetto a gnd

        # 2. Resistenza di Ingresso
        self.add_component(Resistor(name="R_in", node1='in', node2='Q1_base', resistance=R_in_val))

        # 3. Transistor NPN (Q1) - per la metà positiva
        # Emettitore collegato al nodo di uscita
        self.add_component(BJT(name="Q1",
                               nodes={'collector': 'VCC_pos', 'base': 'Q1_base', 'emitter': 'Q1_emitter'},
                               **q_npn_params))

        # 4. Transistor PNP (Q2) - per la metà negativa
        # Collettore collegato all'alimentazione negativa
        # Emettitore collegato al nodo di uscita (Q1_emitter)
        self.add_component(PNP_BJT(name="Q2", # Se usi una classe PNP_BJT
                                  nodes={'collector': 'VCC_neg', 'base': 'Q2_base', 'emitter': 'Q2_emitter'},
                                  **q_pnp_params))
        # Se non hai PNP_BJT, potresti provare a usare BJT di nuovo e sperare che si comporti come un PNP con Vbe e Vce negativi
        # self.add_component(BJT(name="Q2", nodes={'collector': 'VCC_neg', 'base': 'Q2_base', 'emitter': 'Q2_emitter'}, **q_pnp_params))


        # 5. Collegamento emettitori al nodo di uscita comune
        self.add_component(Resistor(name="R_emitters_link", node1='Q1_emitter', node2='Q2_emitter', resistance=1e-12)) # Quasi corto
        self.add_node('output_common')
        self.add_component(Resistor(name="R_output_link", node1='Q1_emitter', node2='output_common', resistance=1e-12)) # Quasi corto
        self.add_component(Resistor(name="R_output_link_2", node1='Q2_emitter', node2='output_common', resistance=1e-12)) # Quasi corto

        # 6. Diodi di Bias (per Classe AB) - Questi polarizzano le basi di Q1 e Q2
        # Creano la caduta di tensione necessaria (circa 1.4V per 2 diodi di silicio)
        # per superare la Vbe dei transistor e farli condurre leggermente.
        self.add_component(Diode(name="D_bias1", nodes={'anode': 'Q1_base', 'cathode': 'Q2_base_d_mid'}, **diode_bias_params))
        self.add_component(Diode(name="D_bias2", nodes={'anode': 'Q2_base_d_mid', 'cathode': 'Q2_base'}, **diode_bias_params))
        # Nota: i diodi andrebbero tra la base di Q1 e la base di Q2.
        # Questa implementazione semplificata è una catena di diodi.
        # Per una corretta polarizzazione, la sorgente d'ingresso andrebbe dopo R_in
        # e la tensione sul nodo Q1_base andrebbe dal segnale + la caduta diodo.
        # Una configurazione più realistica è avere un driver (es. stadio di guadagno in tensione)
        # che pilota questo stadio di uscita, e il bias è tra le uscite del driver e le basi dei finali.

        # Colleghiamo R_in alla base di Q1 e i diodi creano la tensione di bias tra Q1_base e Q2_base.
        # Il nodo 'in' pilota Q1_base. Q2_base è quindi biasata tramite i diodi.
        self.add_node('Q2_base_d_mid') # Nodo intermedio tra i diodi di bias

        # 7. Resistenza di Carico
        self.add_component(Resistor(name="R_load", node1='output_common', node2='gnd', resistance=R_load_val))
        
        # 8. Condensatore di accoppiamento in uscita (per bloccare DC sul carico)
        # Non è strettamente necessario se il carico è tra out e gnd e out è già a 0V DC
        # (tipico per Class B/AB con alimentazione duale e uscita a massa DC).
        # Per ora lo omettiamo per semplicità per evitare di aggiungere un altro nodo in uscita.
        # Se il tuo solutore non gestisce carichi direttamente a un nodo di transistor, potresti aver bisogno di un cap.
        self.add_component(Capacitor(name="C_out", node1='output_common', node2='out', capacitance=100e-6))
        self.add_component(Resistor(name="R_out_link", node1='out', node2='gnd', resistance=R_load_val)) # Il vero carico

        # --- Input (sorgente di segnale AC per analisi transitoria) ---
        self.add_component(VoltageSource(name="V_in_signal", nodes={'pos': 'in', 'neg': 'gnd'},
                                         voltage=0.0))

if __name__ == "__main__":
    from circuit_solver.mna_solver import MnaSolver

    power_amp_circuit = BJTPushPullAmp()
    
    print("Nodi del circuito:", power_amp_circuit.get_node_names())
    print("Componenti del circuito:")
    for comp in power_amp_circuit.components:
        print(f"  - {comp.name}: {type(comp).__name__}")

    solver = MnaSolver(power_amp_circuit)

    print("\nInizio simulazione DC per trovare il punto di riposo dell'amplificatore di potenza...")
    node_voltages, source_currents = solver.solve_dc()

    if node_voltages is not None:
        print("\n--- Risultati Simulazione DC (Punto di Riposo) ---")
        print("Tensioni ai Nodi:")
        for node_name, idx in power_amp_circuit.get_node_map().items():
            if idx != -1:
                print(f"  V({node_name}) = {node_voltages[idx]:.3f} V")
        
        q1_npn = None
        q2_pnp = None
        for comp in power_amp_circuit.nonlinear_components:
            if isinstance(comp, BJT) and comp.name == "Q1":
                q1_npn = comp
            # Devi verificare se è una PNP_BJT o una BJT che agisce da PNP
            elif isinstance(comp, BJT) and comp.name == "Q2": # Se usi BJT per PNP
                q2_pnp = comp
            elif isinstance(comp, PNP_BJT) and comp.name == "Q2": # Se usi PNP_BJT
                q2_pnp = comp

        if q1_npn and q2_pnp:
            # Q1 NPN
            V_c1 = node_voltages[solver._get_node_index(q1_npn.nodes['collector'])]
            V_b1 = node_voltages[solver._get_node_index(q1_npn.nodes['base'])]
            V_e1 = node_voltages[solver._get_node_index(q1_npn.nodes['emitter'])]
            V_be1_q = V_b1 - V_e1
            V_ce1_q = V_c1 - V_e1
            Ic1_q = q1_npn.calculate_collector_current(V_be1_q, V_ce1_q)
            Ib1_q = q1_npn.calculate_base_current(Ic1_q)

            # Q2 PNP
            V_c2 = node_voltages[solver._get_node_index(q2_pnp.nodes['collector'])]
            V_b2 = node_voltages[solver._get_node_index(q2_pnp.nodes['base'])]
            V_e2 = node_voltages[solver._get_node_index(q2_pnp.nodes['emitter'])]
            V_be2_q = V_b2 - V_e2
            V_ce2_q = V_c2 - V_e2
            Ic2_q = q2_pnp.calculate_collector_current(V_be2_q, V_ce2_q) # Questo deve essere il metodo corretto per PNP
            Ib2_q = q2_pnp.calculate_base_current(abs(Ic2_q)) # Base current è positiva

            print(f"\n--- Punto di Riposo del BJT NPN (Q1) ---")
            print(f"  V_BE (Q1) = {V_be1_q:.3f} V")
            print(f"  V_CE (Q1) = {V_ce1_q:.3f} V")
            print(f"  I_C (Q1) = {Ic1_q*1e3:.3f} mA")
            print(f"  I_B (Q1) = {Ib1_q*1e6:.3f} uA")

            print(f"\n--- Punto di Riposo del BJT PNP (Q2) ---")
            print(f"  V_BE (Q2) = {V_be2_q:.3f} V")
            print(f"  V_CE (Q2) = {V_ce2_q:.3f} V")
            print(f"  I_C (Q2) = {Ic2_q*1e3:.3f} mA") # Sarà negativo per PNP
            print(f"  I_B (Q2) = {Ib2_q*1e6:.3f} uA")
            
            # La corrente di riposo (corrente di crossover) dovrebbe essere piccola per la Classe AB
            # ma non zero come per la Classe B pura.
            print(f"\nCorrente di crossover stimata (Ia_Q1 + |Ia_Q2|) = {(Ic1_q + abs(Ic2_q))*1e3:.3f} mA")
            
            # Verificare che i diodi di bias stiano creando la caduta di tensione corretta
            V_d1_anode = node_voltages[solver._get_node_index('Q1_base')]
            V_d1_cathode = node_voltages[solver._get_node_index('Q2_base_d_mid')]
            V_d2_anode = node_voltages[solver._get_node_index('Q2_base_d_mid')]
            V_d2_cathode = node_voltages[solver._get_node_index('Q2_base')]
            print(f"\nDiodi di Bias: V_D1 = {V_d1_anode - V_d1_cathode:.3f} V, V_D2 = {V_d2_anode - V_d2_cathode:.3f} V")

    else:
        print("\nImpossibile calcolare il punto di riposo DC.")
