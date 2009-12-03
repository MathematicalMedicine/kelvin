#include "model.h"

void initializeDefaults ();
void readConfigFile (char *config);
void parseCommandLine (int argc, char *argv[]);
void validateConfig ();
void fillConfigDefaults (ModelRange *modelRange, ModelOptions *modelOptions, ModelType *modelType);
void finishConfig (ModelRange *modelRange, ModelType *modelType);
