#include "model_range.h"
#include "model_options.h"
#include "model_type.h"

/* Structure for the configuration parser dispatch table.
 * 'key' is the configuration directive
 * 'parse' is an int-returning function that takes (char**, int, void*)
 * 'hint' is a pointer that will be passed as the last argument to 'parse'
 */
typedef struct {
  char *key;
  int (*parse) (char **, int, void *);
  void *hint;
} st_dispatch;

extern ModelRange modelRange;
extern ModelOptions modelOptions;
extern ModelType modelType;


int readConfigFile (char *file);
void my_readConfigFile (char *config);
void my_parseCommandLine (int argc, char *argv[]);

/* This really belongs in tools or utils */
void permuteLine (char *line);
