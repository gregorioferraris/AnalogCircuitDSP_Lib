// components/VoltageSource.cpp
#include "VoltageSource.h"
#include <cmath> // Per M_PI, std::sin, std::fmod
#include <iostream> // Per messaggi di errore

// Definisci PI se M_PI non è disponibile in cmath su alcuni compilatori
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Implementazione del costruttore della classe VoltageSource.
 * @param name Nome del componente.
 * @param node_names_str Un vettore di stringhe contenente i nomi dei nodi nell'ordine:
 * [output_node, ground_node].
 * @param type Il tipo di forma d'onda da generare.
 * @param initial_voltage Tensione DC iniziale o ampiezza per le forme d'onda AC.
 * @param frequency Frequenza per le forme d'onda AC (Hz).
 * @param phase_offset Offset di fase per le forme d'onda AC (radianti).
 * @param sample_rate Il sample rate audio in Hz (usato per le sorgenti AC).
 */
VoltageSource::VoltageSource(const std::string& name,
                             const std::vector<std::string>& node_names_str,
                             WaveformType type,
                             double initial_voltage, double frequency,
                             double phase_offset, double sample_rate)
    : Component(name, node_names_str), // Chiama il costruttore della classe base
      type_(type),
      voltage_(initial_voltage),
      frequency_(frequency),
      phase_offset_(phase_offset),
      sample_rate_(sample_rate),
      current_phase_(0.0)
{
    // Una sorgente di tensione indipendente richiede una variabile ausiliaria per la sua corrente di ramo.
    // Questo è gestito dal MnaSolver che assegna gli indici delle variabili ausiliarie.
    // Non è necessario chiamare setNumAuxiliaryVariables() qui, ma il MnaSolver deve sapere
    // quanti ne servono per ogni componente.
    std::cout << "VoltageSource " << name_ << " inizializzato (Tipo: "
              << (type_ == DC ? "DC" : (type_ == SINE ? "SINE" : (type_ == SQUARE ? "SQUARE" : "EXTERNAL")))
              << ", Tensione: " << voltage_ << "V, Freq: " << frequency_ << "Hz)." << std::endl;

    // Assicurati che il sample rate sia valido per le forme d'onda AC
    if (sample_rate_ <= 0) sample_rate_ = 44100.0;
    // Assicurati che la frequenza non sia negativa
    if (frequency_ < 0) frequency_ = 0.0;
}

/**
 * @brief Genera la tensione istantanea della sorgente.
 *
 * Questo metodo genera la tensione appropriata in base al `WaveformType` configurato.
 * - `DC`: Restituisce la `voltage_` costante.
 * - `SINE`: Genera un'onda sinusoidale basata su `voltage_` (ampiezza), `frequency_`,
 * `phase_offset_` e `sample_rate_`. La `current_phase_` viene aggiornata per il prossimo campione.
 * - `SQUARE`: Genera un'onda quadra. L'output è `voltage_` per la prima metà del periodo
 * e `-voltage_` per la seconda metà.
 * - `EXTERNAL`: Semplicemente restituisce l'`input_sample` fornito, agendo come un pass-through.
 *
 * @param time Il tempo attuale della simulazione.
 * @return La tensione generata.
 */
double VoltageSource::getInstantaneousVoltage(double time) {
    double output_voltage = 0.0;

    switch (type_) {
        case DC:
            output_voltage = voltage_;
            break;
        case SINE: {
            // Calcola la fase istantanea in radianti
            // Usiamo il tempo assoluto per la generazione, non un phase accumulator incrementale,
            // per evitare accumulo di errori e per essere coerenti con MNA.
            double phase_radians = (2.0 * M_PI * frequency_ * time) + phase_offset_;
            output_voltage = voltage_ * std::sin(phase_radians);
            break;
        }
        case SQUARE: {
            // Calcola la fase istantanea in radianti per la logica dell'onda quadra
            double phase_radians = (2.0 * M_PI * frequency_ * time) + phase_offset_;
            // Normalizza la fase all'intervallo 0-1 per una logica più semplice dell'onda quadra
            double normalized_phase = std::fmod(phase_radians / (2.0 * M_PI), 1.0);
            if (normalized_phase < 0) normalized_phase += 1.0; // Assicurati che sia positivo

            if (normalized_phase < 0.5) { // Prima metà del periodo
                output_voltage = voltage_;
            } else { // Seconda metà del periodo
                output_voltage = -voltage_;
            }
            break;
        }
        case EXTERNAL:
            // Per il tipo EXTERNAL, la tensione dovrebbe essere fornita dall'esterno
            // (es. da un file audio o un altro componente).
            // Questo metodo non ha un parametro input_sample, quindi non può usarlo.
            // In un contesto MNA, una sorgente EXTERNAL potrebbe richiedere un meccanismo
            // per impostare il suo valore in ogni passo temporale.
            // Per ora, restituirà 0.0, a meno che non venga impostato un valore esplicito.
            // Questo caso d'uso specifico potrebbe richiedere un approccio diverso nell'integrazione MNA.
            output_voltage = 0.0; // Placeholder
            std::cerr << "Attenzione: VoltageSource " << name_ << " di tipo EXTERNAL chiamato senza un valore di input esterno." << std::endl;
            break;
        default:
            output_voltage = 0.0;
            break;
    }

    return output_voltage;
}


/**
 * @brief Applica gli "stamps" della sorgente di tensione alla matrice MNA (A) e al vettore (B).
 *
 * Questo metodo introduce una variabile ausiliaria per la corrente che fluisce attraverso la sorgente di tensione.
 * Gli stamps impongono la relazione di tensione V(node_out) - V(node_ground) = V_source.
 *
 * @param num_total_equations Dimensione totale della matrice MNA.
 * @param dt Passo temporale (non usato per questo componente statico).
 * @param x Vettore della soluzione corrente (non usato per lo stamping).
 * @param prev_solution Vettore della soluzione al passo temporale precedente (non usato).
 * @param time Tempo attuale della simulazione. Usato per le sorgenti AC.
 * @param A Riferimento alla matrice MNA.
 * @param B Riferimento al vettore delle sorgenti (RHS).
 */
void VoltageSource::getStamps(
    int num_total_equations, double dt,
    const std::vector<double>& x,
    const std::vector<double>& prev_solution,
    double time,
    std::vector<std::vector<double>>& A,
    std::vector<double>& B
) {
    // Ottieni gli indici globali per i nodi di uscita e ground.
    // Assumiamo che node_ids_ contenga [output_node, ground_node]
    if (node_ids_.size() != 2) {
        std::cerr << "Errore: VoltageSource " << name_ << " si aspetta 2 ID di nodo, ma ne ha " << node_ids_.size() << std::endl;
        return;
    }

    int idx_out = node_ids_[0];
    int idx_gnd = node_ids_[1]; // Questo dovrebbe essere sempre 0 per il ground

    // Ottieni l'indice per la variabile ausiliaria (corrente della sorgente di tensione).
    // L'MnaSolver assegnerà un indice univoco per ogni sorgente di tensione.
    // Assumiamo che `component_id_` sia l'indice della variabile ausiliaria.
    // Questo è un *workaround* temporaneo. MnaSolver dovrebbe gestire questo.
    int aux_idx = component_id_; // Placeholder: l'MnaSolver dovrà assegnare questo correttamente

    // Controlla la validità degli indici
    if (idx_out < 0 || idx_out >= num_total_equations ||
        idx_gnd < 0 || idx_gnd >= num_total_equations ||
        aux_idx < 0 || aux_idx >= num_total_equations) // aux_idx dovrebbe essere > num_nodes
    {
        std::cerr << "Errore: Indice di nodo o ausiliario non valido per VoltageSource " << name_ << std::endl;
        return;
    }

    // Ottieni la tensione istantanea della sorgente
    double V_source = getInstantaneousVoltage(time);

    // --- Applica gli stamps per l'equazione di vincolo della sorgente di tensione ---
    // Equazione: V(idx_out) - V(idx_gnd) = V_source
    // Questa equazione va nella riga corrispondente alla variabile ausiliaria (corrente).
    A[aux_idx][idx_out] += 1.0;
    A[aux_idx][idx_gnd] -= 1.0;
    B[aux_idx] += V_source; // La tensione della sorgente va nel RHS

    // --- Applica i contributi di corrente alle equazioni dei nodi ---
    // La variabile ausiliaria rappresenta la corrente che esce da idx_gnd ed entra in idx_out.
    A[idx_out][aux_idx] += 1.0;
    A[idx_gnd][aux_idx] -= 1.0;
}

/**
 * @brief Aggiorna lo stato interno della sorgente di tensione.
 *
 * Essendo una sorgente di tensione ideale, non ha uno stato interno da aggiornare
 * nel senso di condensatori o induttori. Questo metodo è fornito per soddisfare
 * l'interfaccia della classe base Component.
 *
 * @param current_solution Il vettore della soluzione corrente (non usato).
 * @param prev_solution Il vettore della soluzione precedente (non usato).
 * @param dt Il passo temporale (non usato).
 */
void VoltageSource::updateState(const std::vector<double>& current_solution,
                         const std::vector<double>& prev_solution,
                         double dt) {
    // Nessuno stato interno da aggiornare per una sorgente di tensione ideale.
    // Questo metodo è intenzionalmente lasciato vuoto.
}

// Setter implementations (dal tuo codice originale)
void VoltageSource::setWaveformType(WaveformType type) {
    type_ = type;
    // Reset phase if changing to an AC type to ensure a consistent start
    if (type_ == SINE || type_ == SQUARE) {
        current_phase_ = 0.0;
    }
}

void VoltageSource::setVoltage(double voltage) {
    voltage_ = voltage;
}

void VoltageSource::setFrequency(double freq) {
    if (freq >= 0) { // Frequency cannot be negative
        frequency_ = freq;
        // Reset phase if frequency changes to avoid jumps in waveform
        current_phase_ = 0.0;
    }
}

void VoltageSource::setPhaseOffset(double offset) {
    phase_offset_ = offset;
}

void VoltageSource::setSampleRate(double rate) {
    if (rate > 0) { // Sample rate must be positive
        sample_rate_ = rate;
        // Reset phase if sample rate changes to avoid jumps in waveform
        current_phase_ = 0.0;
    }
}
