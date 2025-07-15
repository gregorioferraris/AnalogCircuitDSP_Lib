# CMakeLists.txt (del tuo simulatore, aggiornato per i sottocircuiti)
cmake_minimum_required(VERSION 3.10)
project(AnalogCircuitSimulator CXX)

# Imposta la versione standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Aggiungi le directory dei sorgenti
# Le directory sono relative alla root del progetto
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/circuit_solver
    ${CMAKE_CURRENT_SOURCE_DIR}/circuit_solver/subcircuits # NUOVA INCLUSIONE
    ${CMAKE_CURRENT_SOURCE_DIR}/components
    # ${CMAKE_CURRENT_SOURCE_DIR}/utils/physical_quantities # Se hai aggiunto le grandezze fisiche
)

# Definisci i file sorgente
set(SOURCE_FILES
    main.cpp
    circuit_solver/Circuit.cpp
    circuit_solver/MnaSolver.cpp
    circuit_solver/subcircuits/Subcircuit.cpp # NUOVO FILE
    # components/Component.cpp # Se hai un Component.cpp
    components/Resistor.cpp
    components/Capacitor.cpp
    components/Inductor.cpp
    components/VoltageSource.cpp
    components/CurrentSource.cpp
    components/Splitter.cpp
    components/Diode.cpp
    components/SchottkyDiode.cpp
    components/ZenerDiode.cpp
    components/MOSFET.cpp
    components/Triode.cpp
    components/Pentode.cpp
    components/RectifierTube.cpp
    components/LED.cpp
    components/JFET.cpp
    components/BJT.cpp
    components/LDR.cpp
    components/SpeakerDriver.cpp
    components/ClosedBoxCabinet.cpp
    components/BassReflexCabinet.cpp
    # Aggiungi qui tutti gli altri file .cpp dei componenti man mano che li crei
    
    # utils/physical_quantities/PhysicalQuantity.cpp # Se hai aggiunto le grandezze fisiche
    # utils/physical_quantities/Voltage.cpp
    # utils/physical_quantities/Current.cpp
    # utils/physical_quantities/Power.cpp
)

# Crea un eseguibile
add_executable(simulator ${SOURCE_FILES})

# Imposta le opzioni di compilazione
target_compile_options(simulator PRIVATE
    -Wall # Abilita tutti gli avvisi
    -Wextra # Abilita avvisi aggiuntivi
    -g # Abilita le informazioni di debug
    -O2 # Ottimizzazione
)

# Istruzioni per la compilazione:
# 1. Assicurati di avere CMake e un compilatore C++ (es. g++ o clang) installati.
# 2. Nella directory radice del progetto (dove si trova questo CMakeLists.txt), crea una directory 'build':
#    mkdir build
# 3. Entra nella directory 'build':
#    cd build
# 4. Esegui CMake per configurare il progetto:
#    cmake ..
# 5. Compila il progetto:
#    cmake --build .
# 6. L'eseguibile 'simulator' (o 'simulator.exe' su Windows) sarà creato nella directory 'build'.
#    Puoi eseguirlo con: ./simulator
