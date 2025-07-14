# components/speaker_driver.py
import numpy as np
from components.component import Component

class SpeakerDriver(Component):
    def __init__(self, name: str, elec_plus_node: str, elec_minus_node: str, mech_velocity_node: str,
                 Re: float, Le: float, Bl: float, Mms: float, Cms: float, Rms: float, Sd: float):
        """
        Inizializza un modello di driver per altoparlante basato su parametri Thiele-Small.
        Utilizza analogie elettriche per le parti meccaniche.
        Args:
            name (str): Nome univoco dell'istanza (es. "SPKR1").
            elec_plus_node (str): Nodo elettrico positivo della bobina mobile.
            elec_minus_node (str): Nodo elettrico negativo della bobina mobile.
            mech_velocity_node (str): Nodo che rappresenta la velocità meccanica del cono (tensione analogica).
                                      Questo nodo si connetterà al cabinet.
            Re (float): Resistenza DC della bobina mobile (Ohm).
            Le (float): Induttanza della bobina mobile (Henry).
            Bl (float): Fattore di forza (Tesla-metri o N/A).
            Mms (float): Massa meccanica in movimento del cono (kg).
            Cms (float): Cedevolezza meccanica della sospensione (m/N).
            Rms (float): Resistenza meccanica delle perdite (Ns/m).
            Sd (float): Area effettiva del cono (m^2).
        """
        super().__init__(name, elec_plus_node, elec_minus_node, mech_velocity_node)
        self.pin_names = ('elec_plus', 'elec_minus', 'mech_velocity_out') # Nomi dei pin per mappatura

        # Parametri elettrici
        self.Re = Re
        self.Le = Le

        # Parametri meccanici (convertiti in analogie elettriche)
        # Mms (massa) -> Induttanza (L_mech)
        # Cms (cedevolezza) -> Capacità (C_mech)
        # Rms (resistenza) -> Resistenza (R_mech)
        self.L_mech = Mms
        self.C_mech = Cms
        self.R_mech = Rms
        self.Sd = Sd # Area effettiva del cono

        # Fattore di accoppiamento elettro-meccanico
        self.Bl = Bl # Agisce come un trasformatore ideale o gyrator

        # Variabili di stato per Le, L_mech, C_mech (se gestite internamente)
        self.v_Le_prev = 0.0 # Tensione ai capi di Le
        self.i_Le_prev = 0.0 # Corrente attraverso Le
        self.v_Lmech_prev = 0.0 # Tensione ai capi di L_mech (velocità meccanica)
        self.i_Lmech_prev = 0.0 # Corrente attraverso L_mech (forza meccanica)
        self.v_Cmech_prev = 0.0 # Tensione ai capi di C_mech
        self.i_Cmech_prev = 0.0 # Corrente attraverso C_mech

    def get_stamps(self, num_total_equations: int, dt: float, current_solution_guess: np.ndarray, prev_solution: np.ndarray, time: float):
        """
        Restituisce i contributi del driver dell'altoparlante alla matrice MNA (stamp_A) e al vettore RHS (stamp_B).
        Questo metodo implementa un modello a elementi concentrati per l'altoparlante.
        
        Il modello prevede:
        1. Parte Elettrica: Re (Resistor) in serie con Le (Inductor).
        2. Parte Meccanica (analogia elettrica): L_mech (Inductor), C_mech (Capacitor), R_mech (Resistor) in parallelo.
        3. Accoppiamento Bl: Una sorgente dipendente che converte la tensione meccanica (velocità) in corrente elettrica (back-EMF)
           e la corrente elettrica (forza) in tensione meccanica. Questo è il "gyrator".

        Per semplicità, qui implementiamo un modello lineare che può essere gestito direttamente in MNA.
        Assumiamo nodi interni per la bobina e per l'accoppiamento.
        """
        stamp_A = np.zeros((num_total_equations, num_total_equations))
        stamp_B = np.zeros(num_total_equations)

        elec_plus_id = self.node_ids['elec_plus']
        elec_minus_id = self.node_ids['elec_minus']
        mech_velocity_out_id = self.node_ids['mech_velocity_out']
        ground_id = 0 # Assumiamo che il ground sia il nodo 0

        # Nodi interni per il modello dell'altoparlante (devono essere gestiti dalla classe Circuit)
        # Per questo esempio, li useremo come indici virtuali, ma in un sistema reale,
        # Circuit dovrebbe assegnare ID a questi nodi interni.
        # Per ora, si assume che il MnaSolver possa gestire nodi "virtuali" o che questi siano mappati.
        # Questo è un modello semplificato che non crea nodi MNA aggiuntivi, ma integra le equazioni.

        # --- Parte Elettrica (Re in serie con Le) ---
        # Re: Resistore tra elec_plus_id e un nodo interno (non esplicitato qui)
        # Le: Induttore tra nodo interno e elec_minus_id
        # Per MNA, possiamo trattare Re e Le come un'impedenza combinata o come elementi separati.
        # Qui, per get_stamps, li consideriamo come elementi separati che contribuiscono.

        # Contributi di Re (Resistore)
        G_Re = 1.0 / self.Re
        # Assumiamo che Re sia tra elec_plus_id e un nodo virtuale 'elec_mid_node'
        # Per semplicità, consideriamo Re come parte della conduttanza totale tra i nodi elettrici.
        # Questo è un placeholder molto semplificato. Un'implementazione corretta richiederebbe un nodo interno.

        # Contributi di Le (Induttore)
        # Per MNA, l'induttore è modellato con una conduttanza equivalente e una sorgente di corrente equivalente.
        G_eq_Le = dt / (2.0 * self.Le)
        
        # --- Parte Meccanica (L_mech, C_mech, R_mech in parallelo) ---
        # Questi si collegano tra mech_velocity_out_id e ground (riferimento meccanico)
        
        # Contributi di R_mech (Resistore)
        G_Rmech = 1.0 / self.R_mech
        if mech_velocity_out_id != 0: stamp_A[mech_velocity_out_id, mech_velocity_out_id] += G_Rmech

        # Contributi di C_mech (Capacitor)
        G_eq_Cmech = 2.0 * self.C_mech / dt
        if mech_velocity_out_id != 0: stamp_A[mech_velocity_out_id, mech_velocity_out_id] += G_eq_Cmech
        V_Cmech_prev = prev_solution[mech_velocity_out_id] - prev_solution[ground_id]
        i_eq_Cmech = G_eq_Cmech * V_Cmech_prev + self.i_Cmech_prev
        if mech_velocity_out_id != 0: stamp_B[mech_velocity_out_id] -= i_eq_Cmech

        # Contributi di L_mech (Inductor)
        G_eq_Lmech = dt / (2.0 * self.L_mech)
        if mech_velocity_out_id != 0: stamp_A[mech_velocity_out_id, mech_velocity_out_id] += G_eq_Lmech
        V_eq_Lmech = self.i_Lmech_prev * (2.0 * self.L_mech / dt) + self.v_Lmech_prev
        i_eq_Lmech = G_eq_Lmech * V_eq_Lmech
        if mech_velocity_out_id != 0: stamp_B[mech_velocity_out_id] -= i_eq_Lmech

        # --- Accoppiamento Bl (Gyrator) ---
        # Questo è il più complesso da implementare direttamente in get_stamps senza nodi ausiliari.
        # Modella la relazione: Forza = Bl * Corrente_elettrica, Tensione_elettrica_indotta = Bl * Velocità_meccanica
        # Nel dominio della frequenza, questo è un trasformatore ideale.
        # In MNA, richiede sorgenti dipendenti dalla tensione/corrente.
        
        # Per ora, per MNA, useremo una semplificazione dove Bl accoppia direttamente.
        # Questo è un modello lineare, quindi può essere aggiunto direttamente ad A e B.
        # La corrente elettrica (tra elec_plus_id e elec_minus_id) crea una forza (corrente nel dominio meccanico)
        # La velocità meccanica (tensione su mech_velocity_out_id) crea una back-EMF (tensione nel dominio elettrico)

        # Sorgente di tensione dipendente dalla velocità (Back-EMF)
        # V_back_emf = Bl * Velocità_cono
        # Questa tensione è in serie con Re e Le.
        # A_matrix[elec_plus_id, mech_velocity_out_id] -= self.Bl # Esempio di accoppiamento

        # Sorgente di corrente dipendente dalla corrente elettrica (Forza)
        # I_force = Bl * Corrente_elettrica
        # Questa corrente viene iniettata nel nodo di velocità meccanica.
        # A_matrix[mech_velocity_out_id, elec_plus_id] += self.Bl # Esempio di accoppiamento

        # L'implementazione completa del gyrator è complessa per get_stamps.
        # Per il pre-calcolo della risposta in frequenza, useremo le funzioni get_impedance.

        return stamp_A, stamp_B

    def update_state(self, v_Le_curr: float, i_Le_curr: float, v_Lmech_curr: float, i_Lmech_curr: float, v_Cmech_curr: float, i_Cmech_curr: float):
        """
        Aggiorna lo stato interno dei componenti dinamici equivalenti.
        Questo metodo dovrebbe essere chiamato dal MnaSolver per i componenti dinamici.
        """
        self.v_Le_prev = v_Le_curr
        self.i_Le_prev = i_Le_curr
        self.v_Lmech_prev = v_Lmech_curr
        self.i_Lmech_prev = i_Lmech_curr
        self.v_Cmech_prev = v_Cmech_curr
        self.i_Cmech_prev = i_Cmech_curr

    def get_electrical_impedance(self, frequency: float) -> complex:
        """
        Calcola l'impedenza elettrica della bobina mobile (Re + jwLe).
        Args:
            frequency (float): Frequenza in Hz.
        Returns:
            complex: Impedenza elettrica.
        """
        omega = 2 * np.pi * frequency
        return self.Re + 1j * omega * self.Le

    def get_mechanical_impedance_analogy(self, frequency: float) -> complex:
        """
        Calcola l'impedenza meccanica equivalente (in analogia elettrica) del cono e sospensione.
        (R_mech in parallelo con L_mech e C_mech).
        Args:
            frequency (float): Frequenza in Hz.
        Returns:
            complex: Impedenza meccanica in analogia elettrica.
        """
        omega = 2 * np.pi * frequency
        
        # Impedenza di R_mech
        Z_Rmech = self.R_mech
        
        # Impedenza di L_mech (massa)
        Z_Lmech = 1j * omega * self.L_mech
        
        # Impedenza di C_mech (cedevolezza)
        Z_Cmech = 1.0 / (1j * omega * self.C_mech)
        
        # Combinazione in parallelo: 1/Z_total = 1/Z_R + 1/Z_L + 1/Z_C
        # Evita divisione per zero a 0 Hz per induttori/condensatori
        Y_total = 0.0
        if Z_Rmech != 0:
            Y_total += 1.0 / Z_Rmech
        if omega != 0: # Evita divisione per zero per L_mech e C_mech a DC
            if Z_Lmech != 0:
                Y_total += 1.0 / Z_Lmech
            if Z_Cmech != 0:
                Y_total += 1.0 / Z_Cmech
        else: # Comportamento a DC: L_mech è corto, C_mech è aperto
            if Z_Rmech != 0:
                Y_total += 1.0 / Z_Rmech # Solo R_mech contribuisce a DC
        
        return 1.0 / Y_total if Y_total != 0 else np.inf

    def get_velocity_transfer_function(self, frequency: float) -> complex:
        """
        Calcola la funzione di trasferimento dalla tensione elettrica in ingresso alla velocità del cono.
        Questo è il cuore dell'accoppiamento elettro-acustico per il driver.
        Assumiamo che l'ingresso sia una tensione V_in applicata ai terminali elettrici.
        La velocità V_mech_out è la tensione sul nodo meccanico.
        
        V_mech_out / V_in = (Bl / (Re + jwLe)) * (1 / (Z_mech + (Bl^2 / (Re + jwLe))))
        Questo è un modello semplificato che non include la radiazione acustica.
        """
        Z_elec = self.get_electrical_impedance(frequency)
        Z_mech_analogy = self.get_mechanical_impedance_analogy(frequency)
        
        # Impedenza elettrica "riflessa" dalla parte meccanica
        Z_reflected = (self.Bl**2) / Z_mech_analogy if Z_mech_analogy != 0 else np.inf
        
        # Corrente attraverso la bobina mobile: I_elec = V_in / (Z_elec + Z_reflected)
        # Forza sul cono: F_mech = Bl * I_elec
        # Velocità del cono: V_mech = F_mech / Z_mech_analogy
        # V_mech / V_in = (Bl / (Z_elec + Z_reflected)) * (1 / Z_mech_analogy)
        # Questo è il guadagno di transconduttanza (velocità per tensione in ingresso)
        
        # Per ottenere la velocità in m/s per V_in in Volt:
        # V_mech_analogy = (Bl * V_in) / (Z_elec * Z_mech_analogy + Bl^2)
        # In analogia elettrica, V_mech_analogy è la tensione sul nodo meccanico.
        
        denominator = (Z_elec * Z_mech_analogy) + (self.Bl**2)
        if denominator == 0:
            return np.inf # Risonanza infinita
        
        # Questo è il guadagno tra V_elettrica e V_meccanica (velocità)
        # Il termine Bl / Z_elec è la transconduttanza iniziale.
        # Il termine Z_mech_analogy / (Z_mech_analogy + Z_reflected_mech) è la divisione di corrente.
        
        # La funzione di trasferimento dalla tensione elettrica alla velocità del cono è:
        # H(f) = (Bl / (Z_elec + Z_reflected)) * (1 / Z_mech_analogy)
        # Oppure, più direttamente, la corrente elettrica è V_in / (Z_elec + Z_reflected)
        # La forza è Bl * I_elec
        # La velocità è Forza / Z_mech_analogy
        # Quindi, V_velocity / V_in = Bl / ( (Re + jwLe) * Z_mech_analogy + Bl^2 )
        
        # Per un modello di driver, la velocità del cono è proporzionale al flusso volumetrico acustico.
        # La velocità del cono (V_mech_out) è la tensione sul nodo 'mech_velocity_out'.
        # La corrente che entra in questo nodo è la forza.
        
        # La funzione di trasferimento dalla tensione elettrica alla velocità del cono è data da:
        # V_cone_velocity / V_electrical_input = Bl / ( (Re + jwLe) * Z_mechanical_analogy + Bl^2 )
        
        # Z_total_elec = Z_elec + (Bl**2 / Z_mech_analogy) # Impedenza elettrica totale vista dall'ingresso
        
        # Corrente elettrica I_elec = V_in / Z_total_elec
        # Forza F_mech = Bl * I_elec
        # Velocità V_mech = F_mech / Z_mech_analogy
        # Quindi, V_mech / V_in = (Bl / Z_total_elec) / Z_mech_analogy = Bl / (Z_total_elec * Z_mech_analogy)
        # Sostituendo Z_total_elec:
        # V_mech / V_in = Bl / ( (Z_elec + Bl**2 / Z_mech_analogy) * Z_mech_analogy )
        # = Bl / ( Z_elec * Z_mech_analogy + Bl**2 )
        
        # Questo è il guadagno dalla tensione elettrica alla velocità meccanica (tensione analogica)
        return self.Bl / denominator if denominator != 0 else np.inf

