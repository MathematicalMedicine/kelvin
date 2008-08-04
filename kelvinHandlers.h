void quitSignalHandler (int signal);
void usr1SignalHandler (int signal);
void termSignalHandler (int signal);
void intSignalHandler (int signal);
void exitKelvin ();
void setupHandlers ();

extern volatile sig_atomic_t statusRequestSignal;

