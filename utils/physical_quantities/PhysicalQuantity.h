// utils#ifndef PHYSICAL_QUANTITY_H
#define PHYSICAL_QUANTITY_H

#include <string>
#include <iostream>
#include <cmath> // Per operazioni matematiche

class PhysicalQuantity {
protected:
    double value;
    std::string unit;

public:
    // Costruttore
    PhysicalQuantity(double val, const std::string& u) : value(val), unit(u) {}

    // Distruttore virtuale per polimorfismo
    virtual ~PhysicalQuantity() = default;

    // Getter
    double getValue() const { return value; }
    const std::string& getUnit() const { return unit; }

    // Setter (se necessario, ma spesso le quantità sono immutabili dopo la creazione)
    void setValue(double val) { value = val; }
    void setUnit(const std::string& u) { unit = u; }

    // Operatori di base (virtuali per permettere specializzazioni)
    virtual PhysicalQuantity operator+(const PhysicalQuantity& other) const {
        if (unit != other.unit) {
            throw std::invalid_argument("Cannot add quantities with different units.");
        }
        return PhysicalQuantity(value + other.value, unit);
    }

    virtual PhysicalQuantity operator-(const PhysicalQuantity& other) const {
        if (unit != other.unit) {
            throw std::invalid_argument("Cannot subtract quantities with different units.");
        }
        return PhysicalQuantity(value - other.value, unit);
    }

    // Operatore di moltiplicazione (generico, da specializzare nelle classi derivate)
    // Non possiamo restituire un tipo specifico qui senza conoscere il risultato
    // Es: Volt * Ampere = Watt. Quindi, questo sarà un placeholder.
    // Le classi derivate implementeranno operator* specifici.
    // virtual PhysicalQuantity operator*(const PhysicalQuantity& other) const = 0;

    // Operatore di divisione (generico, da specializzare)
    // virtual PhysicalQuantity operator/(const PhysicalQuantity& other) const = 0;

    // Stampa la grandezza fisica
    virtual void print() const {
        std::cout << value << " " << unit;
    }

    // Overload dell'operatore << per una stampa più comoda
    friend std::ostream& operator<<(std::ostream& os, const PhysicalQuantity& pq) {
        pq.print();
        return os;
    }
};

#endif // PHYSICAL_QUANTITY_H
