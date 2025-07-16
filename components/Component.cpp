// components/Component.cpp
#include "Component.h"
#include "../utils/Random.h" // Include la nuova utility Random
#include <stdexcept> // For std::runtime_error

Component::Component(const std::string& name, const std::vector<std::string>& nodeNames)
    : name(name), nodeNamesStr(nodeNames), componentId(-1) {
    // Inizializza nodeIdsVector e nodeIdsMap come vuoti.
    // Saranno popolati da Circuit::addComponent.
}

void Component::setNodeIds(const std::vector<int>& ids) {
    if (ids.size() != nodeNamesStr.size()) {
        throw std::runtime_error("Mismatch between node names count and node IDs count for component " + name);
    }
    nodeIdsVector = ids;
    nodeIdsMap.clear(); // Assicurati che l'altra mappa sia vuota
}

void Component::setNodeIds(const std::map<std::string, int>& ids) {
    if (ids.size() != nodeNamesStr.size()) {
        throw std::runtime_error("Mismatch between node names count and node IDs count for component " + name);
    }
    nodeIdsMap = ids;
    nodeIdsVector.clear(); // Assicurati che l'altra mappa sia vuota
}

void Component::setComponentId(int id) {
    componentId = id;
}

