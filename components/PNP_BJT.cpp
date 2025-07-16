// File: PNP_BJT.cpp
#include "PNP_BJT.h"
#include <algorithm> // Per std::max (non direttamente usato nella versione finale, ma utile per logiche generali dei componenti)

/**
 * @brief Implementazione del costruttore della classe PNP_BJT.
 * @param name Nome del componente.
 * @param input_node Nome del nodo di input.
 * @param output_node Nome del nodo di output.
 * @param ground_node Nome del nodo di massa.
 * @param beta Guadagno di corrente.
 * @param saturation_current Corrente di saturazione inversa.
 * @param thermal_voltage Tensione termica.
 * @param ideality_factor Fattore di idealità.
 */
PNP_BJT::PNP_BJT(const std::string& name, const std::string& input_node, const std::string& output_node,
                 const std::string& ground_node, double beta, double saturation_current,
                 double thermal_voltage, double ideality_factor)
    : Component(name, input_node, output_node, ground_node),
      beta_(beta),
      saturation_current_(saturation_current),
      thermal_voltage_(thermal_voltage),
      ideality_factor_(ideality_factor)
{
    // Assicurati che i parametri siano ragionevoli
    if (beta_ <= 0) beta_ = 1.0;
    if (saturation_current_ <= 0) saturation_current_ = 1e-15; // Valore positivo molto piccolo
    if (thermal_voltage_ <= 0) thermal_voltage_ = 0.001; // Valore positivo piccolo
    if (ideality_factor_ <= 0) ideality_factor_ = 1.0;
}

/**
 * @brief Elabora un singolo campione audio attraverso il modello BJT PNP.
 *
 * Questo modello calcola la corrente di collettore (Ic) basandosi sul campione di input,
 * che è interpretato come la tensione Base-Emettitore (Vbe).
 *
 * Viene utilizzato il modello semplificato di Ebers-Moll per la corrente di collettore (Ic):
 * Ic = Beta * Is * (exp(Vbe / (n * Vt)) - 1)
 *
 * Per un transistor PNP, per la conduzione, Vbe (input_sample) deve essere negativa
 * (cioè, la tensione di Base è inferiore alla tensione di Emettitore).
 *
 * @param input_sample Il campione di input, interpretato come la tensione Base-Emettitore (Vbe).
 * @return La corrente di Collettore (Ic) simulata.
 */
double PNP_BJT::process(double input_sample) {
    // input_sample è Vbe (tensione Base-Emettitore)
    double Vbe = input_sample;

    // Calcola l'argomento per la funzione esponenziale
    double exp_arg = Vbe / (ideality_factor_ * thermal_voltage_);

    // Limita l'argomento esponenziale per prevenire overflow/underflow per valori estremi.
    // Per Vbe grandi e positivi, exp_arg può essere molto grande, portando a overflow.
    // Per Vbe grandi e negativi, exp_arg può essere molto piccolo, portando a underflow (exp(molto_negativo) ~ 0).
    // Il range tipico per exp_arg nei modelli BJT è solitamente entro +/- 20-30 per valori pratici di Vbe.
    const double MAX_EXP_ARG = 20.0; // Corrisponde a Vbe ~ 0.52V per n=1, Vt=0.026V
    const double MIN_EXP_ARG = -20.0; // Corrisponde a Vbe ~ -0.52V per n=1, Vt=0.026V

    if (exp_arg > MAX_EXP_ARG) {
        exp_arg = MAX_EXP_ARG;
    } else if (exp_arg < MIN_EXP_ARG) {
        exp_arg = MIN_EXP_ARG;
    }

    // Calcola il termine del diodo per la giunzione base-emettitore
    double Ib_diode_term = saturation_current_ * (std::exp(exp_arg) - 1.0);

    // Per un PNP, la corrente di collettore è proporzionale alla corrente di base.
    // Stiamo modellando Ic direttamente da Vbe.
    // Un modello semplificato per Ic è Beta * Ib_diode_term.
    // Tuttavia, per i PNP, la conduzione avviene quando Vbe è negativa.
    // Il termine esponenziale (exp(Vbe / (n * Vt)) - 1) sarà negativo quando Vbe è negativo.
    // Quindi, Ic sarà negativo, il che è coerente con la corrente che esce dal collettore.
    // Per l'audio, spesso ci interessa la magnitudine del segnale o il suo effetto sulla tensione.
    // Se l'output deve essere positivo per un effetto di "guadagno", potremmo prendere il valore assoluto o invertire.
    // Per ora, atteniamoci al modello fisico.

    double collector_current = beta_ * Ib_diode_term;

    // Importante: In un circuito reale, la corrente di collettore è anche limitata dal resistore di collettore
    // e dalla tensione di alimentazione. Questo modello calcola solo la corrente ideale basata su Vbe.
    // Per gli effetti audio, questa corrente verrebbe poi convertita di nuovo in una tensione attraverso un resistore di carico.
    // Se l'output deve essere limitato (es. per il clipping), ciò avverrebbe in uno stadio successivo
    // o limitando direttamente questo valore.

    return collector_current;
}
