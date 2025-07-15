// components/Component.cpp
#include "Component.h"
#include <random> // Per std::random_device

// Inizializzazione del generatore di numeri casuali statico per la classe Component.
// Usa std::random_device{}() per un seed più casuale in produzione.
// Per test e riproducibilità, puoi usare un numero fisso (es. 0).
std::mt19937 Component::rng(std::random_device{}()); 
// std::mt19937 Component::rng(0); // Per riproducibilità nei test
