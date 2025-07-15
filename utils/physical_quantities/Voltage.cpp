// utils/physical_quantities/Voltage.cpp
#include "Voltage.h"
#include "Current.h" // Inclusione necessaria per la definizione di Current
#include "Power.h"   // Inclusione necessaria per la definizione di Power

// Definizione dell'operatore di moltiplicazione Voltage * Current = Power
Power Voltage::operator*(const Current& other) const {
    return Power(value * other.value);
}
