# utils/physical_quantities.py

class PhysicalQuantity:
    """
    Rappresenta una grandezza fisica con un valore numerico e un'unità di misura.
    """
    def __init__(self, value: float, unit: str):
        if not isinstance(value, (int, float)):
            raise TypeError("Il valore deve essere un numero.")
        if not isinstance(unit, str) or not unit:
            raise ValueError("L'unità deve essere una stringa non vuota.")
        
        self.value = float(value)
        self.unit = unit.strip() # Rimuovi spazi bianchi

    def __str__(self):
        """Rappresentazione stringa della grandezza fisica (es. "100.0 Ohm")."""
        return f"{self.value} {self.unit}"

    def __repr__(self):
        """Rappresentazione per il debug."""
        return f"PhysicalQuantity(value={self.value}, unit='{self.unit}')"

    # --- Operazioni Aritmetiche Base ---
    # Per semplicità, queste operazioni si concentrano sulla compatibilità di unità
    # e sulla derivazione di unità semplici. Per una gestione completa di unità complesse
    # (es. Ohm * Farad = Secondi), sarebbe necessaria una libreria di gestione unità più sofisticata.

    def __add__(self, other):
        if isinstance(other, PhysicalQuantity):
            if self.unit != other.unit:
                raise ValueError(f"Impossibile sommare grandezze con unità diverse: {self.unit} vs {other.unit}")
            return PhysicalQuantity(self.value + other.value, self.unit)
        return PhysicalQuantity(self.value + other, self.unit) # Assume 'other' è un numero senza unità

    def __sub__(self, other):
        if isinstance(other, PhysicalQuantity):
            if self.unit != other.unit:
                raise ValueError(f"Impossibile sottrarre grandezze con unità diverse: {self.unit} vs {other.unit}")
            return PhysicalQuantity(self.value - other.value, self.unit)
        return PhysicalQuantity(self.value - other, self.unit)

    def __mul__(self, other):
        if isinstance(other, PhysicalQuantity):
            # Esempio di derivazione di unità (molto semplificato)
            new_unit = f"{self.unit}*{other.unit}"
            # Potresti aggiungere logica per semplificare unità comuni (es. V*A=W, Ohm*F=s)
            return PhysicalQuantity(self.value * other.value, new_unit)
        return PhysicalQuantity(self.value * other, self.unit)

    def __truediv__(self, other):
        if isinstance(other, PhysicalQuantity):
            new_unit = f"{self.unit}/{other.unit}"
            return PhysicalQuantity(self.value / other.value, new_unit)
        return PhysicalQuantity(self.value / other, self.unit)

    def __rmul__(self, other): # Per permettere 2 * PhysicalQuantity(...)
        return self.__mul__(other)

    def __rtruediv__(self, other): # Per permettere 1 / PhysicalQuantity(...)
        if isinstance(other, PhysicalQuantity):
            new_unit = f"{other.unit}/{self.unit}"
            return PhysicalQuantity(other.value / self.value, new_unit)
        return PhysicalQuantity(other / self.value, f"1/{self.unit}")

    # --- Comparazioni ---
    def __lt__(self, other):
        if isinstance(other, PhysicalQuantity):
            if self.unit != other.unit:
                raise ValueError(f"Impossibile confrontare grandezze con unità diverse: {self.unit} vs {other.unit}")
            return self.value < other.value
        return self.value < other

    def __le__(self, other):
        if isinstance(other, PhysicalQuantity):
            if self.unit != other.unit:
                raise ValueError(f"Impossibile confrontare grandezze con unità diverse: {self.unit} vs {other.unit}")
            return self.value <= other.value
        return self.value <= other

    def __gt__(self, other):
        if isinstance(other, PhysicalQuantity):
            if self.unit != other.unit:
                raise ValueError(f"Impossibile confrontare grandezze con unità diverse: {self.unit} vs {other.unit}")
            return self.value > other.value
        return self.value > other

    def __ge__(self, other):
        if isinstance(other, PhysicalQuantity):
            if self.unit != other.unit:
                raise ValueError(f"Impossibile confrontare grandezze con unità diverse: {self.unit} vs {other.unit}")
            return self.value >= other.value
        return self.value >= other

    def __eq__(self, other):
        if isinstance(other, PhysicalQuantity):
            return self.value == other.value and self.unit == other.unit
        # Potresti voler confrontare solo il valore se l'altra non ha unità
        return self.value == other # Assume 'other' è un numero senza unità

    def __ne__(self, other):
        return not self.__eq__(other)

    # --- Metodi per l'accesso al valore (se si vuole estrarlo senza unità) ---
    def get_value(self) -> float:
        return self.value

    def get_unit(self) -> str:
        return self.unit
