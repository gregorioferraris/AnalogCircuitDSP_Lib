# CMakeLists.txt (del tuo simulatore, aggiornato per le grandezze fisiche)
cmake_minimum_required(VERSION 3.10)
project(AnalogCircuitSimulator CXX)

# Imposta la versione standard C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Aggiungi le directory dei sorgenti
# Le directory sono relative alla root del progetto
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/circuit_solver
    ${CMAKE_CURRENT_SOURCE_DIR}/components
    ${CMAKE_CURRENT_SOURCE_DIR}/utils/physical_quantities # NUOVA INCLUSIONE
)

# Definisci i file sorgente
set(SOURCE_FILES
    main.cpp
    circuit_solver/Circuit.cpp
    circuit_solver/MnaSolver.cpp
    # components/Component.cpp # Component.h non ha un .cpp separato se le funzioni sono inline o template
    components/Resistor.cpp
    components/Capacitor.cpp
    components/Inductor.cpp
    components/VoltageSource.cpp
    components/CurrentSource.cpp
    components/Splitter.cpp
    # components/Diode.cpp # Se hai già implementato Diode.cpp, includilo
    
    # NUOVI FILE PER LE GRANDEZZE FISICHE
    utils/physical_quantities/PhysicalQuantity.cpp
    utils/physical_quantities/Voltage.cpp
    utils/physical_quantities/Current.cpp
    utils/physical_quantities/Power.cpp
    # Aggiungi qui tutti gli altri file .cpp dei componenti man mano che li crei
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
