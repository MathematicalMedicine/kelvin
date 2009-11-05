#include "model_range.h"
#include "model_options.h"
#include "model_type.h"

extern ModelRange *modelRange;
extern ModelOptions *modelOptions;
extern ModelType *modelType;

void initializeDefaults ();
void readConfigFile (char *config);
void parseCommandLine (int argc, char *argv[]);
void validateConfig ();
void fillConfigDefaults (ModelRange *modelRange, ModelOptions *modelOptions, ModelType *modelType);
void finishConfig (ModelRange *modelRange, ModelType *modelType);
