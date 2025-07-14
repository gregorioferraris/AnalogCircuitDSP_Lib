# utils/speaker_group_utility.py
import numpy as np
from components.speaker_driver import SpeakerDriver # Assicurati che questo import sia corretto

def create_equivalent_speaker_driver(
    base_driver: SpeakerDriver,
    num_drivers: int,
    wiring_config: str,
    name: str = "Equivalent_Speaker_Group"
) -> SpeakerDriver:
    """
    Crea un singolo oggetto SpeakerDriver equivalente che rappresenta un gruppo di driver identici
    connessi in una specifica configurazione.

    Args:
        base_driver (SpeakerDriver): Un'istanza del driver base (es. un singolo cono da 12 pollici).
                                     Si assumono tutti i driver nel gruppo identici a questo.
        num_drivers (int): Il numero totale di driver nel gruppo (es. 4 per un 4x12).
        wiring_config (str): La configurazione di cablaggio ('series', 'parallel', 'series_parallel_4xN').
                             'series_parallel_4xN' è per 4 driver (2 in serie, poi 2 di questi gruppi in parallelo).
        name (str): Nome per il driver equivalente risultante.

    Returns:
        SpeakerDriver: Un nuovo oggetto SpeakerDriver con i parametri Thiele-Small equivalenti.
    """
    if num_drivers <= 0:
        raise ValueError("Il numero di driver deve essere maggiore di zero.")
    if num_drivers == 1:
        return base_driver # Se c'è un solo driver, restituisci quello base

    Re_eq = base_driver.Re
    Le_eq = base_driver.Le
    Bl_eq = base_driver.Bl
    Mms_eq = base_driver.Mms
    Cms_eq = base_driver.Cms
    Rms_eq = base_driver.Rms
    Sd_eq = base_driver.Sd

    if wiring_config == 'series':
        Re_eq = num_drivers * base_driver.Re
        Le_eq = num_drivers * base_driver.Le
        Bl_eq = num_drivers * base_driver.Bl # Bl si somma in serie per la forza totale
        # Mms, Cms, Rms rimangono invariati per singolo cono, ma Sd si somma
        Sd_eq = num_drivers * base_driver.Sd
    elif wiring_config == 'parallel':
        Re_eq = base_driver.Re / num_drivers
        Le_eq = base_driver.Le / num_drivers
        Bl_eq = base_driver.Bl # Bl non cambia per singolo cono, ma la corrente si divide
        Sd_eq = num_drivers * base_driver.Sd
    elif wiring_config == 'series_parallel_4xN':
        if num_drivers != 4:
            raise ValueError("La configurazione 'series_parallel_4xN' è intesa per 4 driver.")
        
        # Tipico 4x12: 2 driver in serie, poi 2 di questi gruppi in parallelo.
        # Calcola prima i parametri per un gruppo di 2 in serie
        Re_series_pair = 2 * base_driver.Re
        Le_series_pair = 2 * base_driver.Le
        Bl_series_pair = 2 * base_driver.Bl
        
        # Poi combina 2 di questi gruppi in parallelo
        Re_eq = Re_series_pair / 2
        Le_eq = Le_series_pair / 2
        Bl_eq = Bl_series_pair # Bl rimane lo stesso per il gruppo in parallelo
        Sd_eq = num_drivers * base_driver.Sd # L'area totale è sempre la somma

    else:
        raise ValueError(f"Configurazione di cablaggio non supportata: {wiring_config}")

    # I parametri meccanici Mms, Cms, Rms si riferiscono al singolo cono
    # e non cambiano direttamente con il cablaggio elettrico.
    # Tuttavia, nel contesto di un "driver equivalente" che rappresenta l'intero gruppo,
    # la massa totale Mms_eq dovrebbe essere la somma delle masse dei singoli driver.
    # E la conformità Cms_eq dovrebbe essere la conformità di un singolo driver divisa per N.
    # Questo è un dettaglio importante per la precisione del modello equivalente.
    # Per un modello equivalente che si muove come un'unica entità:
    Mms_eq = num_drivers * base_driver.Mms
    Cms_eq = base_driver.Cms / num_drivers # La cedevolezza totale del sistema diminuisce
    Rms_eq = num_drivers * base_driver.Rms # La resistenza meccanica totale aumenta

    return SpeakerDriver(
        name=name,
        elec_plus_node="elec_in_plus", # Questi nodi sono solo placeholder per il driver equivalente
        elec_minus_node="elec_in_minus",
        mech_velocity_node="mech_vel",
        Re=Re_eq,
        Le=Le_eq,
        Bl=Bl_eq,
        Mms=Mms_eq,
        Cms=Cms_eq,
        Rms=Rms_eq,
        Sd=Sd_eq
    )

