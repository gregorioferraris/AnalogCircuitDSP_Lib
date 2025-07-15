// utils/physical_quantities/Voltage.h
#ifndef VOLTAGE_H
#define VOLTAGE_H

#include "PhysicalQuantity.h"
#include "Current.h" // Forward declaration o inclusione per operator*
// #include "Resistance.h" // Per operator/

class Current; // Forward declaration

class Voltage : public PhysicalQuantity {
public:
    Voltage(double val) : PhysicalQuantity(val, "V") {}

    // Operatore di moltiplicazione: Voltage * Current = Power
    Power operator*(const Current& other) const;

    // Operatore di divisione: Voltage / Resistance = Current
    // Current operator/(const Resistance& other) const; // Richiede Resistance.h
};

#endif // VOLTAGE_H
