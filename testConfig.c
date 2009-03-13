// gcc -o testConfig testConfig.c config.c -I pedlib -I utils -L lib -lped -lutils -lm

#include "pedlib.h"
ModelOptions modelOptions;
ModelRange modelRange;
ModelType modelType;

#include "sw.h"

int main (int argc, char *argv[] ) {
  readConfigFile (argv[1]);
}
