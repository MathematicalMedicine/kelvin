#include "model_range.h"
#include "model_options.h"
#include "model_type.h"

extern ModelRange modelRange;
extern ModelOptions modelOptions;
extern ModelType modelType;

int readConfigFile (char *file);
void my_readConfigFile (char *config);
void parseCommandLine (int argc, char *argv[]);
void validateConfig ();
