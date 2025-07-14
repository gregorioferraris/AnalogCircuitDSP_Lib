# utils/ngspice_lib_parser.py

import re
import os

class NgspiceLibParser:
    """
    Utility per analizzare i file .lib di NGSPICE ed estrarre i parametri dei modelli dei componenti.
    Supporta modelli comuni come DIODE, NPN, PNP, NMOS, PMOS, JFET, TRIODE.
    """
    def __init__(self):
        self.models = {} # Dizionario per memorizzare i modelli estratti { 'model_name': { 'type': 'DIODE', 'params': { ... } } }

    def parse_lib_file(self, file_path: str):
        """
        Analizza un file .lib di NGSPICE e popola il dizionario dei modelli.
        Args:
            file_path (str): Percorso del file .lib da analizzare.
        Raises:
            FileNotFoundError: Se il file non esiste.
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File .lib non trovato: {file_path}")

        print(f"Analisi del file .lib: {file_path}")
        with open(file_path, 'r') as f:
            content = f.read()

        # Regex per trovare blocchi .MODEL
        # Cattura il nome del modello, il tipo e i parametri
        # Esempio: .MODEL D1N4148 D (Is=10f N=1.0 Rs=16m Cjo=2p M=0.33 Vj=0.7 TT=12n Bv=100 Ibv=10u)
        model_pattern = re.compile(r'\.MODEL\s+([a-zA-Z0-9_]+)\s+([a-zA-Z0-9_]+)\s*\(([^)]*)\)', re.IGNORECASE | re.DOTALL)

        for match in model_pattern.finditer(content):
            model_name = match.group(1).upper() # Nome del modello (es. D1N4148)
            model_type = match.group(2).upper() # Tipo di componente (es. D, NPN, NMOS)
            params_str = match.group(3)        # Stringa dei parametri (es. Is=10f N=1.0 ...)

            params = self._parse_parameters_string(params_str)
            
            self.models[model_name] = {
                'type': model_type,
                'params': params
            }
            print(f"  Trovato modello: {model_name} (Tipo: {model_type})")

    def _parse_parameters_string(self, params_str: str) -> dict:
        """
        Analizza una stringa di parametri (es. "Is=10f N=1.0 Rs=16m") in un dizionario.
        Gestisce i suffissi di scala (p, n, u, m, k, meg, g, t).
        Args:
            params_str (str): Stringa di parametri.
        Returns:
            dict: Dizionario dei parametri { 'param_name': value }.
        """
        parsed_params = {}
        # Regex per trovare coppie chiave=valore, gestendo spazi e virgole
        param_pattern = re.compile(r'([a-zA-Z_][a-zA-Z0-9_]*)\s*=\s*([+\-]?\d*\.?\d+(?:[eE][+\-]?\d+)?(?:[pnumkgtMG]?)?)')
        
        scale_factors = {
            'p': 1e-12, 'n': 1e-9, 'u': 1e-6, 'm': 1e-3,
            'k': 1e3, 'meg': 1e6, 'g': 1e9, 't': 1e12
        }

        for match in param_pattern.finditer(params_str):
            key = match.group(1).strip()
            value_str = match.group(2).strip()
            
            # Estrai il suffisso di scala
            suffix = ''
            if value_str and value_str[-1].lower() in scale_factors:
                suffix = value_str[-1].lower()
                value_str = value_str[:-1]
            elif value_str and value_str[-3:].lower() == 'meg': # Special case for 'meg'
                suffix = 'meg'
                value_str = value_str[:-3]

            try:
                value = float(value_str)
                if suffix:
                    value *= scale_factors[suffix]
                parsed_params[key] = value
            except ValueError:
                print(f"  Warning: Impossibile convertire il valore '{value_str}' per il parametro '{key}'. Ignorato.")
        return parsed_params

    def get_model_params(self, model_name: str) -> dict:
        """
        Restituisce i parametri di un modello specifico.
        Args:
            model_name (str): Il nome del modello (es. "D1N4148").
        Returns:
            dict: Dizionario contenente 'type' e 'params' del modello.
        Raises:
            ValueError: Se il modello non è stato trovato.
        """
        model_info = self.models.get(model_name.upper())
        if model_info is None:
            raise ValueError(f"Modello '{model_name}' non trovato nel file .lib analizzato.")
        return model_info

# --- Esempio di Utilizzo ---
if __name__ == "__main__":
    # --- Passo 1: Crea un file .lib di esempio per il test ---
    # Questo è un esempio semplificato. I file .lib reali sono molto più complessi.
    example_lib_content = """
    * Modelli di diodi
    .MODEL D1N4148 D (Is=2.52n N=1.752 Rs=0.568 Cjo=4p M=0.333 Vj=0.7 TT=12n Bv=100 Ibv=10u)
    .MODEL 1N4001 D (Is=18.8n N=2.000 Rs=0.043 Cjo=39.8p M=0.333 Vj=0.4 TT=4.3u Bv=50 Ibv=10u)

    * Modelli di transistor NPN
    .MODEL 2N2222A NPN (Is=10f Bf=200 Vaf=100 Ikf=0.1 Ne=1.5 Ise=10f Rc=1 Cjc=8p Mjc=0.33 Vjc=0.7 Tf=200p Tr=10n)

    * Modelli di MOSFET
    .MODEL NMOS_TEST NMOS (Vto=0.7 Kp=20u Gamma=0.5 Lambda=0.05)
    .MODEL PMOS_TEST PMOS (Vto=-0.7 Kp=10u Gamma=0.4 Lambda=0.06)

    * Modelli di JFET
    .MODEL JFET_TEST NJF (Vto=-2.0 Beta=1m Lambda=0.01)

    * Modelli di TRIODE (es. per valvole, spesso modelli proprietari o Koren)
    * Nota: I modelli di valvola non sono standard SPICE e la loro sintassi .MODEL può variare molto.
    * Questo è un esempio concettuale che potrebbe non essere direttamente compatibile con tutti i modelli di valvola.
    .MODEL 12AX7_JJ TRIODE (MU=100 KP=1000 EX=1.5 KG1=1.0 VGO=0.0)
    """
    
    lib_file_name = "example_models.lib"
    with open(lib_file_name, "w") as f:
        f.write(example_lib_content)
    print(f"File '{lib_file_name}' creato per il test.")

    # --- Passo 2: Utilizza il parser ---
    parser = NgspiceLibParser()
    try:
        parser.parse_lib_file(lib_file_name)

        # Accedi ai modelli estratti
        print("\n--- Modelli Estratti ---")
        for model_name, model_info in parser.models.items():
            print(f"Modello: {model_name}, Tipo: {model_info['type']}, Parametri: {model_info['params']}")

        # Esempio di recupero di un modello specifico
        d_model = parser.get_model_params("D1N4148")
        print(f"\nParametri per D1N4148: {d_model}")

        npn_model = parser.get_model_params("2N2222A")
        print(f"Parametri per 2N2222A: {npn_model}")

        triode_model = parser.get_model_params("12AX7_JJ")
        print(f"Parametri per 12AX7_JJ: {triode_model}")

        # Esempio di modello non trovato
        try:
            parser.get_model_params("NON_ESISTENTE")
        except ValueError as e:
            print(f"\nErrore previsto: {e}")

    except FileNotFoundError as e:
        print(e)
    finally:
        # Pulisci il file di esempio
        if os.path.exists(lib_file_name):
            os.remove(lib_file_name)
            print(f"File '{lib_file_name}' rimosso.")

