#include <stdio.h>
#include <stdlib.h>
#include "tpl.h"		/* TPL serialization, open source under
				 * Berkeley license */

char            fileName[128];
char           *tplFormat = "f#f#f#f#f#f#";	/* Six fixed vectors of
						 * doubles */
char           *traitFileFormat = "%s_trait.tpl";

/* I'm being passed data that is both self-referential and not contiguous in
 * memory, so I can't serialize the whole thing in one shot. What a pain! */

void 
saveTrait(char *pedigree, double **lDT)
{
  tpl_node       *tn;

  sprintf(fileName, traitFileFormat, pedigree);
  tn = tpl_map(tplFormat, lDT[0], 275, lDT[1], 275, lDT[2], 275, lDT[3], 275, lDT[4], 275, lDT[5], 275);
  tpl_pack(tn, 0);
  tpl_dump(tn, TPL_FILE, fileName);
  tpl_free(tn);
  return;
}

double        **
restoreTrait(char *pedigree)
{
  tpl_node       *tn;
  FILE           *file;
  double        **lDT;

  if ((lDT = (double **) malloc(6 * sizeof(double *))) == NULL) {
    fprintf(stderr, "In restoreTrait, malloc of 1st dimension failed!\n");
    exit(1);
  }
  sprintf(fileName, traitFileFormat, pedigree);
  tn = tpl_map(tplFormat, lDT[0], 275, lDT[1], 275, lDT[2], 275, lDT[3], 275, lDT[4], 275, lDT[5], 275);
  if ((file = fopen(fileName, "r"))) {
    fclose(file);
    tpl_load(tn, TPL_FILE, fileName);
    tpl_unpack(tn, 0);
    return lDT;
  }
  return NULL;
}

#ifdef TEST_MAIN

#include <math.h>

int 
main()
{
  int             i, j;
  double        **lDT;

  if ((lDT = (double **) malloc(6 * sizeof(double *))) == NULL) {
    fprintf(stderr, "malloc of 1st dimension failed!\n");
    exit(1);
  }
  for (i = 0; i < 6; i++) {
    if ((lDT[i] = (double *) malloc(275 * sizeof(double))) == NULL) {
      fprintf(stderr, "malloc of 2nd dimension failed\n");
      exit(1);
    }
    for (j = 0; j < 275; j++)
      lDT[i][j] = (double) i + sqrt(j + 1);
  }
  saveTrait("test", lDT);
  free(lDT);
  lDT = restoreTrait("test");
  for (i = 0; i < 6; i++)
    for (j = 0; j < 275; j++)
      fprintf(stderr, "%G ", lDT[i][j]);
  fprintf(stderr, "\n");
  return 0;
}

#endif
