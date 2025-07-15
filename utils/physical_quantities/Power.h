// utils/physical_quantities/Power.h
#ifndef POWER_H
#define POWER_H

#include "PhysicalQuantity.h"

class Power : public PhysicalQuantity {
public:
    Power(double val) : PhysicalQuantity(val, "W") {}

    // Power non ha operator* o operator/ che restituiscano tipi semplici in questo contesto.
    // Potresti definire operator/ per ottenere Voltage o Current se lo desideri:
    // Voltage operator/(const Current& other) const;
    // Current operator/(const Voltage& other) const;
};

#endif // POWER_H
