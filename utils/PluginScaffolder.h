// utils/PluginScaffolder.h
#ifndef PLUGIN_SCAFFOLDER_H
#define PLUGIN_SCAFFOLDER_H

#include <string>
#include <vector>
#include <filesystem> // Per std::filesystem::create_directories, ecc.
#include <iostream>   // Per std::cout, std::cerr
#include <fstream>    // Per std::ofstream
#include <stdexcept>  // Per std::invalid_argument

// Namespace per convenienza
namespace fs = std::filesystem;

class PluginScaffolder {
public:
    PluginScaffolder(const std::string& baseOutputDir = "generated_plugins");

    void createPluginBoilerplate(const std::string& pluginType,
                                 const std::string& pluginName,
                                 const std::string& authorName = "Your Name",
                                 const std::string& description = "A custom audio plugin");

private:
    std::string baseOutputDir;

    void _createFile(const fs::path& path, const std::string& content, bool executable = false);
    
    void _createLv2Boilerplate(const fs::path& projectDir, const std::string& pluginName,
                               const std::string& authorName, const std::string& description);
    
    void _createVst3Boilerplate(const fs::path& projectDir, const std::string& pluginName,
                                const std::string& authorName, const std::string& description,
                                bool isInstrument);
    
    void _createClapBoilerplate(const fs::path& projectDir, const std::string& pluginName,
                                const std::string& authorName, const std::string& description);

    // Funzione helper per generare UUID (necessaria per VST3 e CLAP)
    std::string _generateUuid();
};

#endif // PLUGIN_SCAFFOLDER_H

