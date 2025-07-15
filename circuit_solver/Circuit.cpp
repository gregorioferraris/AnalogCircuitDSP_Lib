// circuit_solver/Circuit.cpp
#include "Circuit.h"
#include "subcircuits/Subcircuit.h" // Includi la definizione del sottocircuito
#include <stdexcept> // Per std::runtime_error

Circuit::Circuit() : nextNodeId(1) {
    nodes["0"] = 0; // Il nodo '0' è sempre il ground e ha ID 0
}

int Circuit::addNode(const std::string& node_name) {
    if (nodes.find(node_name) == nodes.end()) {
        nodes[node_name] = nextNodeId;
        nextNodeId++;
    }
    return nodes[node_name];
}

void Circuit::addComponent(std::shared_ptr<Component> component) {
    // Assegna un ID univoco al componente
    component->setComponentId(static_cast<int>(components.size()));

    // Assegna gli ID numerici ai nodi del componente
    std::vector<int> node_ids_for_component;
    for (const std::string& node_name : component->getNodeNamesStr()) {
        node_ids_for_component.push_back(addNode(node_name));
    }
    component->setNodeIds(node_ids_for_component);

    // Aggiungi il componente alla lista generale
    components.push_back(component);

    // Categorizza i componenti che aggiungono variabili ausiliarie
    if (std::shared_ptr<VoltageSource> vs = std::dynamic_pointer_cast<VoltageSource>(component)) {
        voltageSources.push_back(vs);
    }
    if (std::shared_ptr<Splitter> splitter = std::dynamic_pointer_cast<Splitter>(component)) {
        splitters.push_back(splitter);
    }

    std::cout << "Aggiunto componente: " << component->getName() << " (" << typeid(*component).name() << ")" << std::endl;
}

void Circuit::addSubcircuitInstance(const Subcircuit& subcircuit_def,
                                   const std::map<std::string, std::string>& instance_connections) {
    std::cout << "Aggiunta istanza sottocircuito: " << subcircuit_def.getName() << std::endl;

    // Mappa i nodi interni del sottocircuito ai nodi globali del circuito principale.
    // I nodi non mappati esplicitamente saranno considerati interni al sottocircuito
    // e verranno rinominati per evitare conflitti.
    std::map<int, int> internal_to_global_node_map; // Mappa: ID nodo interno -> ID nodo globale
    std::map<std::string, int> internal_name_to_global_id_map; // Mappa: nome nodo interno -> ID nodo globale

    const Circuit& internal_circuit_def = subcircuit_def.getInternalCircuit();
    const auto& internal_node_map = internal_circuit_def.getNodeMap();

    // 1. Mappa le porte esterne del sottocircuito ai nodi del circuito principale
    for (const std::string& port_name : subcircuit_def.getPortNames()) {
        auto it_conn = instance_connections.find(port_name);
        if (it_conn == instance_connections.end()) {
            throw std::runtime_error("Errore: Porta '" + port_name + "' del sottocircuito '" + subcircuit_def.getName() + "' non connessa.");
        }
        const std::string& global_node_name = it_conn->second;

        // Se la porta è il ground ("0"), usa l'ID globale 0
        if (port_name == "0") {
            if (global_node_name != "0") {
                 std::cerr << "Attenzione: Porta '0' del sottocircuito '" << subcircuit_def.getName() << "' connessa a '" << global_node_name << "'. Il ground del sottocircuito è sempre ID 0." << std::endl;
            }
            internal_to_global_node_map[0] = 0;
            internal_name_to_global_id_map["0"] = 0;
        } else {
            // Aggiungi o ottieni l'ID del nodo globale
            int global_id = addNode(global_node_name);
            int internal_id = internal_node_map.at(port_name); // Ottieni l'ID interno della porta
            internal_to_global_node_map[internal_id] = global_id;
            internal_name_to_global_id_map[port_name] = global_id;
        }
    }

    // 2. Aggiungi i componenti interni del sottocircuito al circuito principale
    //    Rinomina i nodi interni non connessi per evitare conflitti e assegna nuovi ID globali.
    for (const auto& internal_comp_ptr : internal_circuit_def.getComponents()) {
        std::vector<std::string> new_node_names_str;
        std::vector<int> new_node_ids;

        // Crea una copia del componente per evitare di modificare l'originale nella definizione
        std::shared_ptr<Component> new_comp = std::shared_ptr<Component>(internal_comp_ptr->clone()); // Richiede un metodo clone() in Component

        // Se Component non ha un metodo clone(), dovrai ricreare il componente
        // in base al suo tipo (dynamic_pointer_cast e new ComponentType(...))
        // Per ora, assumiamo un metodo clone() virtuale in Component.
        // Se non hai clone(), dovrai fare un lungo if-else if per ogni tipo di componente:
        /*
        if (std::shared_ptr<Resistor> r = std::dynamic_pointer_cast<Resistor>(internal_comp_ptr)) {
            new_comp = std::make_shared<Resistor>(r->getName(), "n1", "n2", r->resistance); // Nomi nodi placeholder
        } else if (std::shared_ptr<Capacitor> c = std::dynamic_pointer_cast<Capacitor>(internal_comp_ptr)) {
            new_comp = std::make_shared<Capacitor>(c->getName(), "n1", "n2", c->capacitance);
        }
        // ... e così via per tutti i tipi di componenti
        */

        // Rinomina il componente per l'istanza (opzionale ma utile per debug)
        new_comp->name = subcircuit_def.getName() + "_" + internal_comp_ptr->getName();

        for (const std::string& internal_node_name : internal_comp_ptr->getNodeNamesStr()) {
            int internal_node_id = internal_node_map.at(internal_node_name);

            // Se il nodo interno è una porta esterna già mappata
            if (internal_to_global_node_map.count(internal_node_id)) {
                new_node_ids.push_back(internal_to_global_node_map.at(internal_node_id));
                new_node_names_str.push_back(nodes_by_id.at(internal_to_global_node_map.at(internal_node_id))); // Assuming a reverse map or iterator
            } else {
                // Nodo interno non mappato esternamente, crea un nuovo nome unico
                std::string unique_internal_node_name = subcircuit_def.getName() + "_" + internal_comp_ptr->getName() + "_node_" + internal_node_name;
                int global_id = addNode(unique_internal_node_name);
                new_node_ids.push_back(global_id);
                new_node_names_str.push_back(unique_internal_node_name);
                internal_to_global_node_map[internal_node_id] = global_id; // Mappa anche i nodi interni non-port
                internal_name_to_global_id_map[internal_node_name] = global_id;
            }
        }
        new_comp->setNodeIds(new_node_ids);
        // Ora aggiungi il componente (con i nodi globali mappati) al circuito principale
        // Nota: addComponent assegnerà un nuovo componentId e gestirà categorizzazione
        addComponent(new_comp);
    }
}
