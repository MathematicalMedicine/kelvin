#include "model_range.h"
#include "model_options.h"
#include "model_type.h"

extern ModelRange modelRange;
extern ModelOptions modelOptions;
extern ModelType modelType;

int readConfigFile (char *file);
void my_readConfigFile (char *config);
void my_parseCommandLine (int argc, char *argv[]);

/* This really belongs in tools or utils */
void permuteLine (char *line);
