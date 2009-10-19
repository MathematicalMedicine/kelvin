#include <signal.h>

void quitSignalHandler (int ourSignal);
void usr1SignalHandler (int ourSignal);
void termSignalHandler (int ourSignal);
void intSignalHandler (int ourSignal);
void exitKelvin ();
void setupHandlers ();

extern volatile sig_atomic_t statusRequestSignal;

