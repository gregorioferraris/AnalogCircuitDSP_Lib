# components/voltage_source.py

class VoltageSource:
    """
    Rappresenta una sorgente di tensione ideale (DC per ora).
    Aggiunge una riga/colonna aggiuntiva nella matrice MNA per la corrente della sorgente.
    """
    def __init__(self, name: str, nodes: dict, voltage: float = 0.0):
        """
        Inizializza una sorgente di tensione.

        Args:
            name (str): Il nome unico della sorgente di tensione (es. 'V1', 'VCC').
            nodes (dict): Un dizionario con i nomi dei nodi:
                          {'pos': 'nome_nodo_positivo', 'neg': 'nome_nodo_negativo'}.
            voltage (float): La tensione in Volt fornita dalla sorgente (V_pos - V_neg).
                             Per la DC analysis, è un valore costante.
                             Per l'AC/Transient, questo può essere l'ampiezza o il valore DC.
        """
        if not isinstance(name, str) or not name:
            raise ValueError("Il nome della sorgente di tensione deve essere una stringa non vuota.")
        if not isinstance(nodes, dict) or 'pos' not in nodes or 'neg' not in nodes:
            raise ValueError("I nodi della sorgente di tensione devono includere 'pos' e 'neg'.")
        if not isinstance(voltage, (int, float)):
            raise ValueError("La tensione deve essere un numero.")

        self.name = name
        self.nodes = nodes
        self.voltage = float(voltage)
        self.index = -1  # Indice assegnato dal solutore MNA per la variabile di corrente

    def set_nodes(self, pos_node: str, neg_node: str):
        """
        Imposta i nomi dei nodi a cui la sorgente di tensione è collegata.
        """
        self.nodes['pos'] = pos_node
        self.nodes['neg'] = neg_node

    def get_node_names(self) -> dict:
        """
        Restituisce i nomi dei nodi a cui la sorgente di tensione è collegata.
        """
        return self.nodes

    def get_voltage(self) -> float:
        """
        Restituisce il valore della tensione della sorgente.
        """
        return self.voltage
    
    def set_mna_index(self, index: int):
        """
        Imposta l'indice della variabile di corrente associata a questa sorgente di tensione
        nella matrice MNA. Questo indice sarà >= num_nodes.
        """
        self.index = index

    def __str__(self):
        return f"VoltageSource(name='{self.name}', nodes={self.nodes}, voltage={self.voltage:.2f}V)"

# Esempio di utilizzo
if __name__ == "__main__":
    v_source = VoltageSource(name="V_test", nodes={'pos': 'node_A', 'neg': 'gnd'}, voltage=5.0)
    print(v_source)
    print(f"Tensione: {v_source.get_voltage()}V")
    print(f"Nodi: {v_source.get_node_names()}")

    v_source_batt = VoltageSource(name="Battery", nodes={'pos': 'VCC', 'neg': 'node_X'}, voltage=9.0)
    print(v_source_batt)
