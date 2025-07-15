// utils/physical_quantities/Current.cpp
#include "Current.h"
#include "Voltage.h" // Inclusione necessaria per la definizione di Voltage
#include "Power.h"   // Inclusione necessaria per la definizione di Power

// Definizione dell'operatore di moltiplicazione Current * Voltage = Power
Power Current::operator*(const Voltage& other) const {
    return Power(value * other.value);
}
