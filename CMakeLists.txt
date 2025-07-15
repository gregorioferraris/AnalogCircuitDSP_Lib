# CMakeLists.txt (per compilare l'utility PluginScaffolder)
cmake_minimum_required(VERSION 3.10)
project(PluginScaffolderUtility CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Includi le directory dei sorgenti
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/utils
)

# Definisci i file sorgente dell'utility
set(SOURCE_FILES
    main.cpp
    utils/PluginScaffolder.cpp
)

# Crea l'eseguibile dell'utility
add_executable(scaffolder ${SOURCE_FILES})

# Imposta le opzioni di compilazione
target_compile_options(scaffolder PRIVATE
    -Wall
    -Wextra
    -g
    -O2
)

# Istruzioni per la compilazione di questa utility:
# 1. Assicurati di avere CMake e un compilatore C++ (es. g++ o clang) installati.
# 2. Nella directory radice del progetto (dove si trova questo CMakeLists.txt), crea una directory 'build':
#    mkdir build
# 3. Entra nella directory 'build':
#    cd build
# 4. Esegui CMake per configurare il progetto:
#    cmake ..
# 5. Compila il progetto:
#    cmake --build .
# 6. L'eseguibile 'scaffolder' (o 'scaffolder.exe' su Windows) sarà creato nella directory 'build'.
#    Puoi eseguirlo con: ./scaffolder
#    Questo creerà le cartelle dei plugin nella directory 'generated_plugins'.
