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

    // NUOVO Esempio 5: Genera un plugin Audio Unit (AU - solo macOS)
    // Questo genererà i file .h, .mm e Info.plist.
    // La compilazione è possibile SOLO su macOS con Xcode.
    scaffolder.createPluginBoilerplate(
        "au",
        "MacToneShaperAU",
        "MacAudioDev",
        "A tone shaping effect for macOS Audio Units."
    );

    // NUOVO Esempio 6: Genera un plugin LADSPA (semplice, senza GUI)
    scaffolder.createPluginBoilerplate(
        "ladspa",
        "SimpleGainLADSPA",
        "LinuxAudioDev",
        "A basic gain control LADSPA plugin."
    );

    // NUOVO Esempio 7: Genera un plugin AAX (Pro Tools)
    // ATTENZIONE: Questo boilerplate NON COMPILA senza l'SDK AAX proprietario di Avid.
    // È incluso solo per mostrare la struttura concettuale.
    scaffolder.createPluginBoilerplate(
        "aax",
        "ProToolsDistortionAAX",
        "ProAudioDev",
        "A distortion effect for Pro Tools (AAX)."
    );

    std::cout << "\nGenerazione dei boilerplate completata. Controlla la cartella 'generated_plugins'." << std::endl;
    std::cout << "Per ogni plugin generato, segui le istruzioni nel rispettivo 'CMakeLists.txt' o 'Makefile' per compilarlo." << std::endl;
    std::cout << "Successivamente, potrai integrare la tua logica di simulazione MNA C++ nel core di elaborazione del plugin." << std::endl;

    return 0;
}

