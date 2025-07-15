// utils/PluginScaffolder.cpp
#include "PluginScaffolder.h"
#include <sstream> // Per std::stringstream
#include <iomanip> // Per std::hex, std::setw, std::setfill
#include <random>  // Per std::random_device, std::mt19937, std::uniform_int_distribution

// Inizializzazione statica del generatore di numeri casuali per UUID
namespace {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 15); // Per generare cifre esadecimali
    std::uniform_int_distribution<> dis8_11(8, 11); // Per la cifra specifica del formato UUID
}

PluginScaffolder::PluginScaffolder(const std::string& baseOutputDir)
    : baseOutputDir(baseOutputDir) {
    fs::create_directories(baseOutputDir);
    std::cout << "Base directory for plugins: '" << fs::absolute(baseOutputDir) << "'" << std::endl;
}

void PluginScaffolder::_createFile(const fs::path& path, const std::string& content, bool executable) {
    fs::create_directories(path.parent_path());
    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        std::cerr << "Error: Could not create file " << path << std::endl;
        return;
    }
    ofs << content;
    ofs.close();
    if (executable) {
        // Set executable permissions (Unix-like systems)
        fs::perms p = fs::status(path).permissions();
        fs::permissions(path, p | fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec);
    }
    std::cout << "  Created: " << path << std::endl;
}

std::string PluginScaffolder::_generateUuid() {
    std::stringstream ss;
    for (int i = 0; i < 32; ++i) {
        if (i == 8 || i == 12 || i == 16 || i == 20) {
            ss << "-";
        }
        if (i == 12) { // Set the 13th character to '4' for UUID version 4
            ss << "4";
        } else if (i == 16) { // Set the 17th character to '8', '9', 'a', or 'b'
            ss << std::hex << dis8_11(gen);
        } else {
            ss << std::hex << dis(gen);
        }
    }
    return ss.str();
}

void PluginScaffolder::createPluginBoilerplate(const std::string& pluginType,
                                              const std::string& pluginName,
                                              const std::string& authorName,
                                              const std::string& description) {
    std::string typeLower = pluginType;
    std::transform(typeLower.begin(), typeLower.end(), typeLower.begin(), ::tolower);

    std::string sanitizedPluginName = pluginName;
    std::replace(sanitizedPluginName.begin(), sanitizedPluginName.end(), ' ', '_');

    fs::path projectDir = fs::path(baseOutputDir) / sanitizedPluginName;
    fs::create_directories(projectDir);
    std::cout << "\nCreating project '" << pluginName << "' (" << pluginType << ") in '" << projectDir << "'..." << std::endl;

    if (typeLower == "lv2") {
        _createLv2Boilerplate(projectDir, pluginName, authorName, description);
    } else if (typeLower == "vst3") {
        _createVst3Boilerplate(projectDir, pluginName, authorName, description, false);
    } else if (typeLower == "vsti") {
        _createVst3Boilerplate(projectDir, pluginName, authorName, description, true);
    } else if (typeLower == "clap") {
        _createClapBoilerplate(projectDir, pluginName, authorName, description);
    } else {
        throw std::invalid_argument("Unsupported plugin type: " + pluginType + ". Supported: 'lv2', 'vst3', 'vsti', 'clap'.");
    }
    
    std::cout << "Boilerplate for '" << pluginName << "' (" << pluginType << ") generated successfully." << std::endl;
    std::cout << "Next steps: Compile the C++ code and integrate your processing logic." << std::endl;
}

void PluginScaffolder::_createLv2Boilerplate(const fs::path& projectDir, const std::string& pluginName,
                                            const std::string& authorName, const std::string& description) {
    std::string sanitizedPluginName = pluginName;
    std::replace(sanitizedPluginName.begin(), sanitizedPluginName.end(), ' ', '_');

    fs::path bundleDir = projectDir / (sanitizedPluginName + ".lv2");
    fs::create_directories(bundleDir);
    fs::path srcDir = bundleDir / "src";
    fs::create_directories(srcDir);

    std::string pluginUri = "http://" + authorName.lower().replace(' ', '') + ".com/" + pluginName.lower().replace(' ', '');

    // manifest.ttl
    std::string manifestContent = R"(@prefix lv2: <http://lv2plug.in/ns/lv2core#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .

<)" + pluginUri + R"(> a lv2:Plugin ;
    lv2:binary <)" + sanitizedPluginName + R"(.so> ;
    rdfs:seeAlso <)" + sanitizedPluginName + R"(.ttl> .
)";
    _createFile(bundleDir / "manifest.ttl", manifestContent);

    // plugin.ttl
    std::string pluginTtlContent = R"(@prefix lv2: <http://lv2plug.in/ns/lv2core#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix doap: <http://usefulinc.com/ns/doap#> .
@prefix pset: <http://lv2plug.in/ns/ext/presets#> .

<)" + pluginUri + R"(> a lv2:Plugin ;
    lv2:requiredFeature lv2:isLive ;
    lv2:optionalFeature lv2:hardRTCapable ;
    lv2:port [
        a lv2:AudioPort ,
          lv2:InputPort ;
        lv2:index 0 ;
        lv2:symbol "input" ;
        lv2:name "Input"
    ] , [
        a lv2:AudioPort ,
          lv2:OutputPort ;
        lv2:index 1 ;
        lv2:symbol "output" ;
        lv2:name "Output"
    ] , [
        a lv2:ControlPort ,
          lv2:InputPort ;
        lv2:index 2 ;
        lv2:symbol "gain" ;
        lv2:name "Gain" ;
        lv2:default 0.0 ;
        lv2:minimum -20.0 ;
        lv2:maximum 20.0 ;
        lv2:portProperty lv2:logarithmic ;
    ] ;
    doap:name ")" + pluginName + R"(" ;
    doap:maintainer [
        doap:name ")" + authorName + R"(" ;
        doap:homepage <http://example.com> ;
        doap:mbox <mailto:your.email@example.com>
    ] ;
    doap:description ")" + description + R"(" ;
    lv2:minorVersion 0 ;
    lv2:microVersion 1 .

# Example preset
<urn:)" + pluginName.lower().replace(' ', '_') + R"(:preset:default>
    a pset:Preset ;
    lv2:appliesTo <)" + pluginUri + R"(> ;
    rdfs:label "Default Preset" ;
    lv2:port [
        lv2:symbol "gain" ;
        pset:value 0.0
    ] .
)";
    _createFile(bundleDir / (sanitizedPluginName + ".ttl"), pluginTtlContent);

    // plugin.h
    std::string pluginHContent = R"(#ifndef )" + sanitizedPluginName.upper() + R"(_H
#define )" + sanitizedPluginName.upper() + R"(_H

#include <lv2/lv2plug.in/ns/lv2core/lv2.h>
#include <cmath> // For log, exp, etc.

// Plugin URI (must match manifest.ttl and plugin.ttl)
#define )" + sanitizedPluginName.upper() + R"(_URI ")" + pluginUri + R"("

class )" + sanitizedPluginName + R"_Plugin {
public:
    )" + sanitizedPluginName + R"_Plugin(double sampleRate);
    ~)" + sanitizedPluginName + R"_Plugin();

    // Host callbacks
    void activate();
    void deactivate();
    void process(const float* input, float* output, unsigned long n_samples);
    void setControl(LV2_URID control_urid, float value); // For atom ports or custom controls

private:
    double sampleRate;
    float* input_port;
    float* output_port;
    float* gain_port; // Example control port

    // --- START: Area for your core processing logic ---
    // Here you can translate your Python logic (e.g., from PluginProcessor, ExternalEffectWrapper) to C++.
    // Example:
    // float current_gain_value;
    // float delay_buffer[44100 * 2]; // Example: 2-second delay buffer
    // int delay_write_idx;

    // If you calculated FIR/IIR coefficients in Python, you would use them here:
    // float fir_coefficients[NUM_TAPS];
    // float fir_history[NUM_TAPS - 1]; // Filter state

    // Example processing function
    void applyProcessing(const float* input, float* output, unsigned long n_samples);
    // --- END: Area for your core processing logic ---
};

#endif // )" + sanitizedPluginName.upper() + R"(_H
)";
    _createFile(srcDir / (sanitizedPluginName + ".h"), pluginHContent);

    // plugin.cpp
    std::string pluginCppContent = R"(#include ")" + sanitizedPluginName + R"(.h"
#include <lv2/lv2plug.in/ns/lv2core/lv2.h>
#include <lv2/lv2plug.in/ns/ext/urid/urid.h> // For LV2_URID_Map, if needed
#include <lv2/lv2plug.in/ns/ext/atom/atom.h> // For Atom, if needed
#include <iostream> // For debug

// LV2 callback functions
static LV2_Handle
instantiate(const LV2_Descriptor* descriptor,
            double                    sample_rate,
            LV2_URID_Map* map,
            LV2_Feature* const* features)
{
    )" + sanitizedPluginName + R"_Plugin* plugin = new )" + sanitizedPluginName + R"_Plugin(sample_rate);
    // Initialize URID mapper if needed
    return plugin;
}

static void
cleanup(LV2_Handle instance)
{
    delete static_cast<)" + sanitizedPluginName + R"_Plugin*>(instance);
}

static void
connect_port(LV2_Handle instance,
             uint32_t     port,
             void* data)
{
    )" + sanitizedPluginName + R"_Plugin* plugin = static_cast<)" + sanitizedPluginName + R"_Plugin*>(instance);
    switch (port) {
    case 0: // Input
        plugin->input_port = static_cast<float*>(data);
        break;
    case 1: // Output
        plugin->output_port = static_cast<float*>(data);
        break;
    case 2: // Gain
        plugin->gain_port = static_cast<float*>(data);
        break;
    default:
        break;
    }
}

static void
activate(LV2_Handle instance)
{
    static_cast<)" + sanitizedPluginName + R"_Plugin*>(instance)->activate();
}

static void
deactivate(LV2_Handle instance)
{
    static_cast<)" + sanitizedPluginName + R"_Plugin*>(instance)->deactivate();
}

static void
run(LV2_Handle instance, uint32_t n_samples)
{
    )" + sanitizedPluginName + R"_Plugin* plugin = static_cast<)" + sanitizedPluginName + R"_Plugin*>(instance);
    plugin->process(plugin->input_port, plugin->output_port, n_samples);
}

// Plugin Constructor
)" + sanitizedPluginName + R"_Plugin::)" + sanitizedPluginName + R"_Plugin(double sr) : sampleRate(sr)
{
    input_port = nullptr;
    output_port = nullptr;
    gain_port = nullptr;
    // Initialize your processing logic here
    // current_gain_value = 0.0;
    // delay_write_idx = 0;
    // std::fill(delay_buffer, delay_buffer + sizeof(delay_buffer)/sizeof(delay_buffer[0]), 0.0f);
}

// Plugin Destructor
)" + sanitizedPluginName + R"_Plugin::~)" + sanitizedPluginName + R"_Plugin()
{
    // Clean up resources here
}

void )" + sanitizedPluginName + R"_Plugin::activate()
{
    // Reset plugin state on activation
    // Example:
    // current_gain_value = 0.0;
    // delay_write_idx = 0;
    // std::fill(delay_buffer, delay_buffer + sizeof(delay_buffer)/sizeof(delay_buffer[0]), 0.0f);
}

void )" + sanitizedPluginName + R"_Plugin::deactivate()
{
    // Clean up plugin state on deactivation
}

void )" + sanitizedPluginName + R"_Plugin::process(const float* input, float* output, unsigned long n_samples)
{
    // Retrieve gain parameter value (if connected)
    // float current_gain = (gain_port) ? *gain_port : 0.0f;
    // current_gain_value = current_gain; // Update internal state if necessary

    // --- START: Area for your core processing logic ---
    // This is the function that will be executed for each audio block.
    // Here you can call your applyProcessing function or insert the logic directly.
    applyProcessing(input, output, n_samples);
    // --- END: Area for your core processing logic ---
}

// Example implementation of the processing logic
void )" + sanitizedPluginName + R"_Plugin::applyProcessing(const float* input, float* output, unsigned long n_samples)
{
    float current_gain_val = (gain_port) ? std::pow(10.0f, *gain_port / 20.0f) : 1.0f; // Convert dB to linear

    for (unsigned long i = 0; i < n_samples; ++i) {
        output[i] = input[i] * current_gain_val; // Simple gain
        // Here you would insert your MNA/DSP logic translated from Python.
        // For example, applying an FIR filter:
        // float filtered_sample = 0.0f;
        // for (int k = 0; k < NUM_TAPS; ++k) {
        //     filtered_sample += fir_coefficients[k] * get_history_sample(i - k); // Pseudo-code
        // }
        // output[i] = filtered_sample;
    }
}

// Method to update control parameters (if there are non-audio parameters)
void )" + sanitizedPluginName + R"_Plugin::setControl(LV2_URID control_urid, float value)
{
    // Handle control parameters here, if they are not directly connected to a port
    // (e.g., parameters that don't change for every block)
}

// The LV2 plugin "factory"
static const LV2_Descriptor descriptor = {
    )" + sanitizedPluginName.upper() + R"(_URI,
    instantiate,
    connect_port,
    run,
    deactivate,
    activate,
    cleanup,
    NULL // extension_data
};

// LV2 entry point function
LV2_SYMBOL_EXPORT
const LV2_Descriptor*
lv2_descriptor(uint32_t index)
{
    switch (index) {
    case 0:
        return &descriptor;
    default:
        return NULL;
    }
}
)";
    _createFile(srcDir / (sanitizedPluginName + ".cpp"), pluginCppContent);

    // Makefile
    std::string makefileContent = R"(# Makefile for the LV2 plugin
# To compile: make
# To install (into your LV2 system): make install
# To clean: make clean

LV2_CFLAGS ?= $(shell pkg-config --cflags lv2)
LV2_LIBS ?= $(shell pkg-config --libs lv2)

CXXFLAGS += $(LV2_CFLAGS) -fPIC -Wall -Wextra -g -O2
LDFLAGS += $(LV2_LIBS) -shared -Wl,-Bsymbolic

TARGET = )" + sanitizedPluginName + R"(.so
SOURCES = )" + sanitizedPluginName + R"(.cpp
OBJECTS = $(SOURCES:.cpp=.o)

BUNDLE_DIR = )" + sanitizedPluginName + R"(.lv2
INSTALL_DIR = $(HOME)/.lv2/$(BUNDLE_DIR)

all: $(BUNDLE_DIR)/$(TARGET)

$(BUNDLE_DIR)/$(TARGET): $(OBJECTS)
    $(CXX) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJECTS): $(SOURCES)
    $(CXX) $(CXXFLAGS) -c $< -o $@

install: all
    mkdir -p $(INSTALL_DIR)
    cp $(BUNDLE_DIR)/*.ttl $(INSTALL_DIR)
    cp $(BUNDLE_DIR)/$(TARGET) $(INSTALL_DIR)

clean:
    rm -f $(OBJECTS) $(BUNDLE_DIR)/$(TARGET)
)";
    _createFile(bundleDir / "Makefile", makefileContent);
}

void PluginScaffolder::_createVst3Boilerplate(const fs::path& projectDir, const std::string& pluginName,
                                             const std::string& authorName, const std::string& description,
                                             bool isInstrument) {
    std::string sanitizedPluginName = pluginName;
    std::replace(sanitizedPluginName.begin(), sanitizedPluginName.end(), ' ', '_');

    fs::path srcDir = projectDir / "src";
    fs::create_directories(srcDir);

    std::string componentUuid = _generateUuid();
    std::string controllerUuid = _generateUuid();

    // plugin.h
    std::string pluginHContent = R"(#pragma once
#include "public.sdk/source/vst/vstparameters.h"
#include "public.sdk/source/vst/vstaudioeffect.h"
#include "pluginterfaces/vst/ivstaudioprocessor.h"
#include "pluginterfaces/vst/ivstevents.h" // For VSTi

// Component UUID (randomly generated, must be unique)
static const Steinberg::FUID )" + sanitizedPluginName + R"__Component_UID (0x)" + componentUuid.substr(0, 8) + R"(, 0x)" + componentUuid.substr(8, 8) + R"(, 0x)" + componentUuid.substr(16, 8) + R"(, 0x)" + componentUuid.substr(24, 8) + R"();
// Controller UUID (randomly generated, must be unique)
static const Steinberg::FUID )" + sanitizedPluginName + R"__Controller_UID (0x)" + controllerUuid.substr(0, 8) + R"(, 0x)" + controllerUuid.substr(8, 8) + R"(, 0x)" + controllerUuid.substr(16, 8) + R"(, 0x)" + controllerUuid.substr(24, 8) + R"();

namespace )" + sanitizedPluginName + R"Plugin {

class )" + sanitizedPluginName + R"Processor : public Steinberg::Vst::AudioEffect
{
public:
    )" + sanitizedPluginName + R"Processor();
    ~)" + sanitizedPluginName + R"Processor() SMTG_OVERRIDE;

    // --- VST3 Overrides ---
    tresult PLUGIN_API initialize(FUnknown* context) SMTG_OVERRIDE;
    tresult PLUGIN_API terminate() SMTG_OVERRIDE;
    tresult PLUGIN_API setActive(TBool state) SMTG_OVERRIDE;
    tresult PLUGIN_API process(Vst::ProcessData& data) SMTG_OVERRIDE;
    tresult PLUGIN_API setupProcessing(Vst::ProcessSetup& setup) SMTG_OVERRIDE;
    tresult PLUGIN_API setBusArrangements(Vst::SpeakerArrangement* inputs, int32 numInputs, Vst::SpeakerArrangement* outputs, int32 numOutputs) SMTG_OVERRIDE;

    // --- START: Area for your core processing logic ---
    // Here you can translate your Python logic (e.g., from PluginProcessor, ExternalEffectWrapper) to C++.
    // Example:
    // float current_gain_value = 0.0f;
    // float delay_buffer[44100 * 2]; // Example: 2-second delay buffer
    // int delay_write_idx = 0;

    // If you calculated FIR/IIR coefficients in Python, you would use them here:
    // float fir_coefficients[NUM_TAPS];
    // float fir_history[NUM_TAPS - 1]; // Filter state

    // Internal processing method
    void applyProcessing(const float* input_buffer_L, const float* input_buffer_R, float* output_buffer_L, float* output_buffer_R, int num_samples);
    // --- END: Area for your core processing logic ---

protected:
    // Plugin parameters
    Steinberg::Vst::Parameter* gainParameter; // Example parameter
    float gainValue;
};

// --- VSTi Specific ---
#if )" + (isInstrument ? "1" : "0") + R"(
class )" + sanitizedPluginName + R"InstrumentProcessor : public )" + sanitizedPluginName + R"Processor
{
public:
    tresult PLUGIN_API process(Vst::ProcessData& data) SMTG_OVERRIDE;
protected:
    void processEvents(Vst::IEventList* eventList);
};
#endif

} // namespace )" + sanitizedPluginName + R"Plugin
)";
    _createFile(srcDir / (sanitizedPluginName + ".h"), pluginHContent);

    // plugin.cpp
    std::string pluginCppContent = R"(#include ")" + sanitizedPluginName + R"(.h"
#include "public.sdk/source/vst/vstsingleaudioeffect.h"
#include "pluginterfaces/base/ibstream.h"
#include <cmath> // For log, exp, etc.
#include <iostream> // For debug

namespace )" + sanitizedPluginName + R"Plugin {

// --- )" + sanitizedPluginName + R"Processor Implementation ---
)" + sanitizedPluginName + R"Processor::)" + sanitizedPluginName + R"Processor()
: AudioEffect()
{
    // Register parameters
    parameters.addParameter(gainParameter = new Vst::Parameter(USTRING("Gain"), "Gain", "dB", 0.0, -20.0, 20.0, 0.0, Vst::ParameterInfo::kCanAutomate));
    // Initialize your processing logic here
    // std::fill(delay_buffer, delay_buffer + sizeof(delay_buffer)/sizeof(delay_buffer[0]), 0.0f);
    // delay_write_idx = 0;
}

)" + sanitizedPluginName + R"Processor::~)" + sanitizedPluginName + R"Processor()
{
    // Clean up resources here
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Processor::initialize(FUnknown* context)
{
    tresult result = AudioEffect::initialize(context);
    if (result == kResultTrue)
    {
        // Add audio buses
        addAudioInput(USTRING("Stereo In"), Vst::SpeakerArr::kStereo);
        addAudioOutput(USTRING("Stereo Out"), Vst::SpeakerArr::kStereo);
        #if )" + (isInstrument ? "1" : "0") + R"(
        // For instruments, add MIDI input
        addEventInput(USTRING("MIDI In"), 1);
        #endif
    }
    return result;
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Processor::terminate()
{
    return AudioEffect::terminate();
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Processor::setActive(TBool state)
{
    // Called when the plugin is activated/deactivated
    if (state)
    {
        // Reset plugin state
        // std::fill(delay_buffer, delay_buffer + sizeof(delay_buffer)/sizeof(delay_buffer[0]), 0.0f);
        // delay_write_idx = 0;
    }
    return AudioEffect::setActive(state);
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Processor::setupProcessing(Vst::ProcessSetup& setup)
{
    // Configure the processor (e.g., sample rate, max block size)
    return AudioEffect::setupProcessing(setup);
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Processor::setBusArrangements(Vst::SpeakerArrangement* inputs, int32 numInputs, Vst::SpeakerArrangement* outputs, int32 numOutputs)
{
    // Set audio bus arrangements
    return AudioEffect::setBusArrangements(inputs, numInputs, outputs, numOutputs);
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Processor::process(Vst::ProcessData& data)
{
    // Read parameters
    if (data.inputParameterChanges)
    {
        int32 numParamsChanged = data.inputParameterChanges->getParameterCount();
        for (int32 i = 0; i < numParamsChanged; ++i)
        {
            Vst::IParamValueQueue* paramQueue = data.inputParameterChanges->getParameterQueue(i);
            if (paramQueue)
            {
                Vst::ParamValue value;
                int32 sampleOffset;
                int32 numPoints = paramQueue->getPointCount();
                if (numPoints > 0 && paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue)
                {
                    if (paramQueue->getParameterId() == gainParameter->id)
                    {
                        gainValue = static_cast<float>(value);
                    }
                }
            }
        }
    }

    // --- START: Area for your core processing logic ---
    // Process audio blocks
    if (data.inputs->sampleBuffers && data.outputs->sampleBuffers)
    {
        float* input_buffer_L = data.inputs[0].channelBuffers[0];
        float* input_buffer_R = data.inputs[0].channelBuffers[1];
        float* output_buffer_L = data.outputs[0].channelBuffers[0];
        float* output_buffer_R = data.outputs[0].channelBuffers[1];
        int num_samples = data.numSamples;

        applyProcessing(input_buffer_L, input_buffer_R, output_buffer_L, output_buffer_R, num_samples);
    }
    // --- END: Area for your core processing logic ---

    return kResultTrue;
}

void )" + sanitizedPluginName + R"Processor::applyProcessing(const float* input_buffer_L, const float* input_buffer_R, float* output_buffer_L, float* output_buffer_R, int num_samples)
{
    float current_gain_linear = std::pow(10.0f, gainValue / 20.0f); // Convert dB to linear

    for (int i = 0; i < num_samples; ++i)
    {
        output_buffer_L[i] = input_buffer_L[i] * current_gain_linear;
        output_buffer_R[i] = input_buffer_R[i] * current_gain_linear;
        // Here you would insert your MNA/DSP logic translated from Python.
        // For example, applying an FIR filter:
        // output_buffer_L[i] = filter_sample(input_buffer_L[i]); // Call to a filter function
    }
}

// --- VSTi Specific Implementation ---
#if )" + (isInstrument ? "1" : "0") + R"(
tresult PLUGIN_API )" + sanitizedPluginName + R"InstrumentProcessor::process(Vst::ProcessData& data)
{
    // Handle MIDI events
    if (data.inputEvents)
    {
        processEvents(data.inputEvents);
    }
    // Call the base audio processor
    return )" + sanitizedPluginName + R"Processor::process(data);
}

void )" + sanitizedPluginName + R"InstrumentProcessor::processEvents(Vst::IEventList* eventList)
{
    // Implement MIDI event handling logic here
    // Example:
    // int32 numEvents = eventList->getEventCount();
    // for (int32 i = 0; i < numEvents; ++i)
    // {
    //     Vst::Event* event = eventList->getEvent(i);
    //     if (event->type == Vst::Event::kNoteOnEvent)
    //     {
    //         // Handle Note On
    //     }
    //     else if (event->type == Vst::Event::kNoteOffEvent)
    //     {
    //         // Handle Note Off
    //     }
    // }
}
#endif

} // namespace )" + sanitizedPluginName + R"Plugin
)";
    _createFile(srcDir / (sanitizedPluginName + ".cpp"), pluginCppContent);

    // controller.h
    std::string controllerHContent = R"(#pragma once
#include "public.sdk/source/vst/vsteditcontroller.h"

namespace )" + sanitizedPluginName + R"Plugin {

class )" + sanitizedPluginName + R"Controller : public Steinberg::Vst::EditController
{
public:
    tresult PLUGIN_API initialize(FUnknown* context) SMTG_OVERRIDE;
    tresult PLUGIN_API terminate() SMTG_OVERRIDE;
    tresult PLUGIN_API setComponentState(Steinberg::IBStream* state) SMTG_OVERRIDE;
    tresult PLUGIN_API setState(Steinberg::IBStream* state) SMTG_OVERRIDE;
    tresult PLUGIN_API getState(Steinberg::IBStream* state) SMTG_OVERRIDE;
    tresult PLUGIN_API setParamNormalized(Vst::ParamID tag, Vst::ParamValue value) SMTG_OVERRIDE;
    Vst::ParamValue PLUGIN_API getParamNormalized(Vst::ParamID tag) SMTG_OVERRIDE;
protected:
    float gainValue;
};

} // namespace )" + sanitizedPluginName + R"Plugin
)";
    _createFile(srcDir / (sanitizedPluginName + "Controller.h"), controllerHContent);

    // controller.cpp
    std::string controllerCppContent = R"(#include ")" + sanitizedPluginName + R"Controller.h"
#include ")" + sanitizedPluginName + R"(.h" // For UUIDs
#include "pluginterfaces/base/ibstream.h"
#include "public.sdk/source/vst/vstparameters.h"
#include <iostream> // For debug

namespace )" + sanitizedPluginName + R"Plugin {

// --- )" + sanitizedPluginName + R"Controller Implementation ---
tresult PLUGIN_API )" + sanitizedPluginName + R"Controller::initialize(FUnknown* context)
{
    tresult result = EditController::initialize(context);
    if (result == kResultTrue)
    {
        // Add parameters to the controller
        parameters.addParameter(new Vst::Parameter(USTRING("Gain"), "Gain", "dB", 0.0, -20.0, 20.0, 0.0, Vst::ParameterInfo::kCanAutomate));
    }
    return result;
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Controller::terminate()
{
    return EditController::terminate();
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Controller::setComponentState(Steinberg::IBStream* state)
{
    // Load component state (e.g., parameter values)
    if (!state) return kResultFalse;
    // Example:
    // float savedGain = 0.0f;
    // state->read(&savedGain, sizeof(float));
    // setParamNormalized(gainParameter->id, savedGain);
    return kResultTrue;
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Controller::setState(Steinberg::IBStream* state)
{
    // Load controller state (e.g., parameter values)
    if (!state) return kResultFalse;
    // Example:
    // float savedGain = 0.0f;
    // state->read(&savedGain, sizeof(float));
    // setParamNormalized(gainParameter->id, savedGain);
    return kResultTrue;
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Controller::getState(Steinberg::IBStream* state)
{
    // Save controller state
    if (!state) return kResultFalse;
    // Example:
    // state->write(&gainValue, sizeof(float));
    return kResultTrue;
}

tresult PLUGIN_API )" + sanitizedPluginName + R"Controller::setParamNormalized(Vst::ParamID tag, Vst::ParamValue value)
{
    // Update parameters when they are changed by the host or GUI
    tresult result = EditController::setParamNormalized(tag, value);
    if (result == kResultTrue)
    {
        if (tag == parameters.getParameterByTitle(USTRING("Gain"))->id)
        {
            gainValue = static_cast<float>(value);
        }
    }
    return result;
}

Vst::ParamValue PLUGIN_API )" + sanitizedPluginName + R"Controller::getParamNormalized(Vst::ParamID tag)
{
    // Return the normalized value of a parameter
    if (tag == parameters.getParameterByTitle(USTRING("Gain"))->id)
    {
        return gainValue;
    }
    return EditController::getParamNormalized(tag);
}

// The VST3 "factory"
extern "C"
{
    // Function to create the plugin component
    __attribute__((visibility("default"))) Steinberg::IPluginFactory* createInstance (void* /*hostContext*/)
    {
        return new Steinberg::Vst::SingleAudioEffect<
            #if )" + (isInstrument ? "1" : "0") + R"(
                )" + sanitizedPluginName + R"InstrumentProcessor
            #else
                )" + sanitizedPluginName + R"Processor
            #endif
            , )" + sanitizedPluginName + R"Controller
        >(
            )" + sanitizedPluginName + R"__Component_UID,
            )" + sanitizedPluginName + R"__Controller_UID
        );
    }
}

} // namespace )" + sanitizedPluginName + R"Plugin
)";
    _createFile(srcDir / (sanitizedPluginName + "Controller.cpp"), controllerCppContent);

    // CMakeLists.txt (simplified)
    std::string cmakeContent = R"(cmake_minimum_required(VERSION 3.10)
project()" + sanitizedPluginName + R"Plugin CXX)

# Path to VST3 SDK (you must set this!)
# set(VST3_SDK_PATH "/path/to/VST_SDK/VST3_SDK")
# if(NOT EXISTS ${VST3_SDK_PATH})
#    message(FATAL_ERROR "VST3 SDK not found. Download it from Steinberg and update VST3_SDK_PATH.")
# endif()

# Include SDK directories
# include_directories(
#     ${VST3_SDK_PATH}/public.sdk/source
#     ${VST3_SDK_PATH}/public.sdk/source/vst
#     ${VST3_SDK_PATH}/pluginterfaces/vst
#     ${VST3_SDK_PATH}/pluginterfaces/base
#     ${VST3_SDK_PATH}/pluginterfaces/gui
# )

# Define sources
set(SOURCES
    src/)" + sanitizedPluginName + R"(.cpp
    src/)" + sanitizedPluginName + R"Controller.cpp
)

# Create the shared plugin library
add_library()" + sanitizedPluginName + R"( SHARED ${SOURCES})

# Set compile features
target_compile_features()" + sanitizedPluginName + R"( PUBLIC cxx_17)
target_compile_options()" + sanitizedPluginName + R"( PRIVATE -fPIC -Wall -Wextra -g -O2)

# Set output name for VST3 plugin
set_target_properties()" + sanitizedPluginName + R"( PROPERTIES
    SUFFIX ".vst3"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/VST3"
)

# Link necessary libraries (may vary depending on system)
# target_link_libraries()" + sanitizedPluginName + R"( PUBLIC
#     # Add specific VST3 libraries here, e.g., vstgui.lib, base.lib
# )

# Compilation instructions:
# 1. Download the VST3 SDK from Steinberg.
# 2. Update VST3_SDK_PATH in this file.
# 3. Create a 'build' directory in the project root: mkdir build
# 4. Enter the 'build' directory: cd build
# 5. Run CMake: cmake ..
# 6. Compile: cmake --build .
# The .vst3 plugin will be in the 'VST3' subdirectory within the 'build' directory.
)";
    _createFile(projectDir / "CMakeLists.txt", cmakeContent);
}

void PluginScaffolder::_createClapBoilerplate(const fs::path& projectDir, const std::string& pluginName,
                                             const std::string& authorName, const std::string& description) {
    std::string sanitizedPluginName = pluginName;
    std::replace(sanitizedPluginName.begin(), sanitizedPluginName.end(), ' ', '_');

    fs::path srcDir = projectDir / "src";
    fs::create_directories(srcDir);

    std::string clapId = authorName.lower().replace(' ', '') + "." + pluginName.lower().replace(' ', '');

    // plugin.h
    std::string pluginHContent = R"(#pragma once
#include <clap/clap.h>
#include <cmath> // For log, exp, etc.

// Plugin ID (must be unique)
static const clap_plugin_id )" + sanitizedPluginName + R"__CLAP_ID = ")" + clapId + R"(" ;

class )" + sanitizedPluginName + R"Plugin {
public:
    )" + sanitizedPluginName + R"Plugin(const clap_plugin_descriptor* desc, const clap_host* host);
    ~)" + sanitizedPluginName + R"Plugin();

    // CLAP methods
    bool activate(double sample_rate, uint32_t min_frames, uint32_t max_frames);
    void deactivate();
    void process(const clap_process_t* process);
    void on_param_value(const clap_event_param_value_t* param_value);

    // Static functions for CLAP interface
    static const clap_plugin_descriptor* get_descriptor();
    static bool init(const clap_plugin_descriptor* descriptor);
    static void deinit(const clap_plugin_descriptor* descriptor);
    static clap_plugin* create(const clap_plugin_descriptor* descriptor, const clap_host* host);
    static void destroy(clap_plugin* plugin);

private:
    const clap_plugin_descriptor* descriptor;
    const clap_host* host;
    double sampleRate;

    // Parameters
    float gainValue = 0.0f; // Example parameter

    // --- START: Area for your core processing logic ---
    // Here you can translate your Python logic (e.g., from PluginProcessor, ExternalEffectWrapper) to C++.
    // Example:
    // float delay_buffer[44100 * 2]; // Example: 2-second delay buffer
    // int delay_write_idx = 0;

    // Internal processing method
    void applyProcessing(const float* input_L, const float* input_R, float* output_L, float* output_R, int num_samples);
    // --- END: Area for your core processing logic ---
};
)";
    _createFile(srcDir / (sanitizedPluginName + ".h"), pluginHContent);

    // plugin.cpp
    std::string pluginCppContent = R"(#include ")" + sanitizedPluginName + R"(.h"
#include <clap/ext/audio-ports.h>
#include <clap/ext/params.h>
#include <clap/ext/note-ports.h> // For instruments
#include <iostream> // For debug
#include <cstring> // For strcmp

// --- CLAP Plugin Implementation ---

)" + sanitizedPluginName + R"Plugin::)" + sanitizedPluginName + R"Plugin(const clap_plugin_descriptor* desc, const clap_host* h)
: descriptor(desc), host(h)
{
    // Initialize your processing logic here
    // std::fill(delay_buffer, delay_buffer + sizeof(delay_buffer)/sizeof(delay_buffer[0]), 0.0f);
    // delay_write_idx = 0;
}

)" + sanitizedPluginName + R"Plugin::~)" + sanitizedPluginName + R"Plugin()
{
    // Clean up resources here
}

bool )" + sanitizedPluginName + R"Plugin::activate(double sr, uint32_t min_frames, uint32_t max_frames)
{
    sampleRate = sr;
    // Reset plugin state on activation
    // std::fill(delay_buffer, delay_buffer + sizeof(delay_buffer)/sizeof(delay_buffer[0]), 0.0f);
    // delay_write_idx = 0;
    return true;
}

void )" + sanitizedPluginName + R"Plugin::deactivate()
{
    // Clean up plugin state on deactivation
}

void )" + sanitizedPluginName + R"Plugin::process(const clap_process_t* process)
{
    // Read parameters from the event queue
    const clap_input_events_t* in_events = process->in_events;
    uint32_t num_events = in_events->count(in_events);
    for (uint32_t i = 0; i < num_events; ++i)
    {
        const clap_event_header_t* hdr = in_events->get(in_events, i);
        if (hdr->type == CLAP_EVENT_PARAM_VALUE)
        {
            const clap_event_param_value_t* param_value = (const clap_event_param_value_t*)hdr;
            on_param_value(param_value);
        }
        // Handle MIDI events here if it's an instrument
        // if (hdr->type == CLAP_EVENT_NOTE_ON) { ... }
    }

    // --- START: Area for your core processing logic ---
    // Process audio blocks
    uint32_t num_samples = process->frames_count;
    const float* input_L = process->audio_inputs[0].data32[0];
    const float* input_R = process->audio_inputs[0].data32[1];
    float* output_L = process->audio_outputs[0].data32[0];
    float* output_R = process->audio_outputs[0].data32[1];

    applyProcessing(input_L, input_R, output_L, output_R, num_samples);
    // --- END: Area for your core processing logic ---
}

void )" + sanitizedPluginName + R"Plugin::on_param_value(const clap_event_param_value_t* param_value)
{
    if (param_value->param_id == 0) // ID 0 for Gain (example)
    {
        gainValue = param_value->value;
    }
}

void )" + sanitizedPluginName + R"Plugin::applyProcessing(const float* input_L, const float* input_R, float* output_L, float* output_R, int num_samples)
{
    float current_gain_linear = std::pow(10.0f, gainValue / 20.0f); // Convert dB to linear

    for (int i = 0; i < num_samples; ++i)
    {
        output_L[i] = input_L[i] * current_gain_linear;
        output_R[i] = input_R[i] * current_gain_linear;
        // Here you would insert your MNA/DSP logic translated from Python.
    }
}

// --- CLAP Descriptor and Factory ---

static const clap_plugin_descriptor )" + sanitizedPluginName + R"__descriptor = {
    .clap_version = CLAP_VERSION_INIT,
    .id = )" + sanitizedPluginName + R"__CLAP_ID,
    .name = ")" + pluginName + R"(",
    .vendor = ")" + authorName + R"(",
    .url = "http://example.com",
    .manual_url = "",
    .support_url = "",
    .version = "0.1.0",
    .description = ")" + description + R"(",
    .features = (const char*[]){
        CLAP_PLUGIN_FEATURE_AUDIO_EFFECT,
        CLAP_PLUGIN_FEATURE_STEREO,
        // CLAP_PLUGIN_FEATURE_INSTRUMENT, // For instruments
        NULL
    },
};

bool )" + sanitizedPluginName + R"Plugin::init(const clap_plugin_descriptor* descriptor) { return true; }
void )" + sanitizedPluginName + R"Plugin::deinit(const clap_plugin_descriptor* descriptor) {}

clap_plugin* )" + sanitizedPluginName + R"Plugin::create(const clap_plugin_descriptor* descriptor, const clap_host* host)
{
    clap_plugin* plugin = (clap_plugin*)malloc(sizeof(clap_plugin));
    plugin->desc = descriptor;
    plugin->plugin_data = new )" + sanitizedPluginName + R"Plugin(descriptor, host);

    // Basic CLAP plugin callbacks
    plugin->init = [](const clap_plugin* p) { return static_cast<)" + sanitizedPluginName + R"Plugin*>(p->plugin_data)->init(p->desc); };
    plugin->activate = [](const clap_plugin* p, double sr, uint32_t min_f, uint32_t max_f) { return static_cast<)" + sanitizedPluginName + R"Plugin*>(p->plugin_data)->activate(sr, min_f, max_f); };
    plugin->deactivate = [](const clap_plugin* p) { static_cast<)" + sanitizedPluginName + R"Plugin*>(p->plugin_data)->deactivate(); };
    plugin->process = [](const clap_plugin* p, const clap_process_t* process) { static_cast<)" + sanitizedPluginName + R"Plugin*>(p->plugin_data)->process(process); };
    plugin->destroy = [](const clap_plugin* p) { delete static_cast<)" + sanitizedPluginName + R"Plugin*>(p->plugin_data); free(p); };

    // Implement necessary extensions (audio-ports, params, etc.)
    // Example: Audio Ports
    static const clap_plugin_audio_ports_t audio_ports_ext = {
        .count = [](const clap_plugin* p, bool is_input) { return 1; },
        .get = [](const clap_plugin* p, uint32_t index, bool is_input, clap_audio_port_info_t* info) {
            if (index != 0) return false;
            info->id = 0;
            info->name = is_input ? "Input" : "Output";
            info->flags = CLAP_AUDIO_PORT_IS_MAIN;
            info->channel_count = 2; // Stereo
            info->port_type = CLAP_AUDIO_PORT_STEREO;
            return true;
        }
    };
    plugin->get_extension = [](const clap_plugin* p, const char* id) -> const void* {
        if (std::strcmp(id, CLAP_PLUGIN_AUDIO_PORTS) == 0) return &audio_ports_ext;
        // Add other extensions (params, state, gui, etc.)
        return nullptr;
    };

    return plugin;
}

// The CLAP "factory"
static const clap_plugin_factory )" + sanitizedPluginName + R"__factory = {
    .get_descriptor = [](uint32_t index) -> const clap_plugin_descriptor* {
        if (index == 0) return &)" + sanitizedPluginName + R"__descriptor;
        return nullptr;
    },
    .create_plugin = )" + sanitizedPluginName + R"Plugin::create,
};

// CLAP entry point function
CLAP_EXPORT const clap_plugin_factory* clap_plugin_factory()
{
    return &)" + sanitizedPluginName + R"__factory;
}
)";
    _createFile(srcDir / (sanitizedPluginName + ".cpp"), pluginCppContent);

    // CMakeLists.txt (simplified)
    std::string cmakeContent = R"(cmake_minimum_required(VERSION 3.10)
project()" + sanitizedPluginName + R"Plugin CXX)

# Path to CLAP SDK (you must set this!)
# set(CLAP_SDK_PATH "/path/to/clap")
# if(NOT EXISTS ${CLAP_SDK_PATH})
#    message(FATAL_ERROR "CLAP SDK not found. Download it from https://github.com/free-audio/clap and update CLAP_SDK_PATH.")
# endif()

# Include SDK directories
# include_directories(${CLAP_SDK_PATH}/include)

# Define sources
set(SOURCES
    src/)" + sanitizedPluginName + R"(.cpp
)

# Create the shared plugin library
add_library()" + sanitizedPluginName + R"( SHARED ${SOURCES})

# Set compile features
target_compile_features()" + sanitizedPluginName + R"( PUBLIC cxx_17)
target_compile_options()" + sanitizedPluginName + R"( PRIVATE -fPIC -Wall -Wextra -g -O2)

# Set output name for CLAP plugin
set_target_properties()" + sanitizedPluginName + R"( PROPERTIES
    SUFFIX ".clap"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/CLAP"
)

# Link necessary libraries (may vary depending on system)
# target_link_libraries()" + sanitizedPluginName + R"( PUBLIC
#     # Add specific CLAP libraries here
# )

# Compilation instructions:
# 1. Download the CLAP SDK.
# 2. Update CLAP_SDK_PATH in this file.
# 3. Create a 'build' directory in the project root: mkdir build
# 4. Enter the 'build' directory: cd build
# 5. Run CMake: cmake ..
# 6. Compile: cmake --build .
# The .clap plugin will be in the 'CLAP' subdirectory within the 'build' directory.
)";
    _createFile(projectDir / "CMakeLists.txt", cmakeContent);
}

