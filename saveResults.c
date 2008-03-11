#include <stdio.h>
#include "tpl.h"		/* TPL serialization, open source under Berkeley license */

char fileName[128];
char *tplFormat = "f#";
char *traitFileFormat = "%s_trait.tpl";

#define LDT_LENGTH 1650

void saveTrait(char *pedigree, double  **lDT) {
  tpl_node *tn;
  sprintf(fileName, traitFileFormat, pedigree);
  tn = tpl_map(tplFormat, lDT, LDT_LENGTH);
  tpl_pack(tn, 0);
  tpl_dump(tn, TPL_FILE, fileName);
  tpl_free(tn);
  return;
}

double ** restoreTrait(char *pedigree) {
  tpl_node *tn;
  FILE *file;
  double **lDT;
  sprintf(fileName, traitFileFormat, pedigree);
  tn = tpl_map(tplFormat, lDT, LDT_LENGTH);
  if (file = fopen(fileName, "r")) {
    fclose(file);
    tpl_load(tn, TPL_FILE, fileName);
    tpl_unpack(tn, 0);
    return lDT;
  }
  return NULL;
}

#ifdef TEST_MAIN

#include <stdlib.h>
#include <math.h>

int main() {
  int i,j;
  double **lDT;

  if ((lDT = (double **) malloc(6*275*sizeof(double))) == NULL) {
    fprintf(stderr, "malloc failed!\n");
    exit(1);
  }
  for (i=0; i<6; i++)
    for (j=0; j<275; j++)
      lDT[i][j] = (double) i+sqrt(j+1);
  saveTrait("test", lDT);
  free(lDT);
  lDT = restoreTrait("test");
  for (i=0; i<6; i++)
    for (j=0; j<275; j++)
      fprintf(stderr, "%G ", lDT[i][j]);
  fprintf(stderr, "\n");
  return 0;
}

#endif
