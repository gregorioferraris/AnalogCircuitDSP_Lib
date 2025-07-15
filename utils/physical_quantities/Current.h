// utils/physical_quantities/Current.h
#ifndef CURRENT_H
#define CURRENT_H

#include "PhysicalQuantity.h"
#include "Voltage.h" // Forward declaration o inclusione per operator*

class Voltage; // Forward declaration

class Current : public PhysicalQuantity {
public:
    Current(double val) : PhysicalQuantity(val, "A") {}

    // Operatore di moltiplicazione: Current * Voltage = Power
    Power operator*(const Voltage& other) const;

    // Operatore di moltiplicazione: Current * Resistance = Voltage
    // Voltage operator*(const Resistance& other) const; // Richiede Resistance.h
};

#endif // CURRENT_H
