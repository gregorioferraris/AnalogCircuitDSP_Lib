// main.cpp
#include <iostream>
#include "utils/PluginScaffolder.h"

int main() {
    std::cout << "--- Generazione Boilerplate Plugin Audio ---" << std::endl;

    PluginScaffolder scaffolder;

    // Esempio 1: Genera un plugin LV2
    scaffolder.createPluginBoilerplate(
        "lv2",
        "MyVintageDistortion",
        "AudioDev",
        "A vintage-style distortion plugin."
    );

    // Esempio 2: Genera un plugin VST3 (effetto)
    scaffolder.createPluginBoilerplate(
        "vst3",
        "AnalogChorusFX",
        "AudioDev",
        "An analog-modeled chorus effect."
    );

    // Esempio 3: Genera un plugin VSTi (strumento)
    scaffolder.createPluginBoilerplate(
        "vsti",
        "SimpleSynthVSTi",
        "AudioDev",
        "A basic polyphonic synthesizer instrument."
    );

    // Esempio 4: Genera un plugin CLAP
    scaffolder.createPluginBoilerplate(
        "clap",
        "ModulationDelayCLAP",
        "AudioDev",
        "A modern delay with modulation capabilities."
    );

    std::cout << "\nGenerazione dei boilerplate completata. Controlla la cartella 'generated_plugins'." << std::endl;
    std::cout << "Per ogni plugin generato, segui le istruzioni nel rispettivo 'CMakeLists.txt' o 'Makefile' per compilarlo." << std::endl;
    std::cout << "Successivamente, potrai integrare la tua logica di simulazione MNA C++ nel core di elaborazione del plugin." << std::endl;

    return 0;
}

