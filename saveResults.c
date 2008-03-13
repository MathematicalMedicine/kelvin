#include <stdio.h>
#include <stdlib.h>
#include "kelvin.h"
#include "saveResults.h"
#include "tpl.h"		/* TPL serialization, open source under
				 * Berkeley license */

char fileName[128];

/* I'm being passed data that is both self-referential and not contiguous in
 * memory, so I can't serialize the whole thing in one shot. What a pain! */

void
dump_lDT (double **lDT)
{
  int i, j;

  for (i = 0; i < 6; i++)
    for (j = 0; j < 275; j++)
      fprintf (stderr, "%G ", lDT[i][j]);
  fprintf (stderr, "\n");
  return;
}

char *traitTPLFormat = "sf#f#f#f#f#f#";	/* String and six fixed vectors of doubles */
char *traitFileFormat = "%s%s_trait.tpl";

int
saveTrait (char *pedigree, double **lDT)
{
  tpl_node *tn;

  //  fprintf(stderr, "in saveTrait for pedigree %s\n", pedigree);
  sprintf (fileName, traitFileFormat, resultsprefix, pedigree);
  tn =
    tpl_map (traitTPLFormat, &pedigree, lDT[0], 275, lDT[1], 275, lDT[2], 275,
	     lDT[3], 275, lDT[4], 275, lDT[5], 275);
  tpl_pack (tn, 0);
  tpl_dump (tn, TPL_FILE, fileName);
  tpl_free (tn);
  return 0;
}

int
restoreTrait (char *pedigree, double **lDT)
{
  tpl_node *tn;
  FILE *file;
  char *checkPedigree;

  //  fprintf(stderr, "in restoreTrait for pedigree %s\n", pedigree);
  sprintf (fileName, traitFileFormat, resultsprefix, pedigree);
  tn =
    tpl_map (traitTPLFormat, &checkPedigree, lDT[0], 275, lDT[1], 275, lDT[2],
	     275, lDT[3], 275, lDT[4], 275, lDT[5], 275);
  if ((file = fopen (fileName, "r"))) {
    fclose (file);
    tpl_load (tn, TPL_FILE, fileName);
    tpl_unpack (tn, 0);
    if (strcmp (pedigree, checkPedigree) == 0) {
      tpl_free (tn);
      //      fprintf(stderr, "restoreTrait successfully restored trait data\n");
      return 1;
    } else {
      fprintf (stderr,
	       "restoreTrait check of requested pedigree of %s vs restored of %s failed, exiting!\n",
	       pedigree, checkPedigree);
      exit (1);
    }
  }
  //  fprintf(stderr, "restoreTrait found no trait data\n");
  tpl_free (tn);
  return 0;
}

char *markerTPLFormat = "siiA(s)f";	/* String, two ints, array of string and a single float */
char *markerFileFormat = "%s%s_%i_";

int
saveMarker (char *pedigree, int chromosome, int markerCount,
	    char **markerNames, double *mDT)
{
  tpl_node *tn;
  int i;
  char *markerName;

  //  fprintf(stderr, "in saveMarker for pedigree %s, chromosome %d w/%d markers of value %G\n",
  //      pedigree, chromosome, markerCount, *mDT);
  sprintf (fileName, markerFileFormat, resultsprefix, pedigree, chromosome);
  for (i = 0; i < markerCount; i++) {
    strcat (fileName, markerNames[i]);
    strcat (fileName, "_");
  }
  strcat (fileName, "marker.tpl");
  tn =
    tpl_map (markerTPLFormat, &pedigree, &chromosome, &markerCount,
	     &markerName, mDT);
  tpl_pack (tn, 0);

  for (i = 0; i < markerCount; i++) {
    markerName = markerNames[i];
    tpl_pack (tn, 1);
  }

  tpl_dump (tn, TPL_FILE, fileName);
  tpl_free (tn);

  return 0;
}

int
restoreMarker (char *pedigree, int chromosome, int markerCount,
	       char **markerNames, double *mDT)
{
  int i, checkChromosome, checkMarkerCount;
  tpl_node *tn;
  FILE *file;
  char *checkPedigree;
  char *markerName;

  //  fprintf(stderr, "in restoreMarker for pedigree %s, chromosome %d w/%d markers\n", 
  //      pedigree, chromosome, markerCount);
  sprintf (fileName, markerFileFormat, resultsprefix, pedigree, chromosome);
  for (i = 0; i < markerCount; i++) {
    strcat (fileName, markerNames[i]);
    strcat (fileName, "_");
  }
  strcat (fileName, "marker.tpl");
  tn =
    tpl_map (markerTPLFormat, &checkPedigree, &checkChromosome,
	     &checkMarkerCount, &markerName, mDT);
  if ((file = fopen (fileName, "r"))) {
    fclose (file);
    tpl_load (tn, TPL_FILE, fileName);
    tpl_unpack (tn, 0);
    i = 0;
    markerName = markerNames[i];
    while (tpl_unpack (tn, 1) > 0)
      markerName = markerNames[i++];

    if ((strcmp (pedigree, checkPedigree) == 0)
	&& (checkChromosome == chromosome)
	&& (checkMarkerCount == markerCount)) {
      tpl_free (tn);
      //      fprintf(stderr, "restoreMarker successfully restored marker data of value %G\n", *mDT);
      return 1;
    } else {
      fprintf (stderr,
	       "restoreMarker check of pedigree/marker/count of %s/%d/%d vs %s/%d/%d failed, exiting!\n",
	       pedigree, chromosome, markerCount, checkPedigree,
	       checkChromosome, checkMarkerCount);
      exit (1);
    }
  }
  //  fprintf(stderr, "restoreMarker found no marker data\n");
  tpl_free (tn);
  return 0;
}

char *alternativeTPLFormat = "siff#f#f#f#f#f#";	/* String and six fixed vectors of doubles */
char *alternativeFileFormat = "%s%s_%i_%G_alternative.tpl";

int
saveAlternative (char *pedigree, int chromosome, double traitPosition,
		 double **lDT)
{
  tpl_node *tn;

  //  fprintf(stderr, "in saveAlternative for pedigree %s, chromosome %d, trait position %d\n",
  //      pedigree, chromosome, (int) traitPosition);
  sprintf (fileName, alternativeFileFormat, resultsprefix, pedigree,
	   chromosome, traitPosition);
  tn =
    tpl_map (alternativeTPLFormat, &pedigree, &chromosome, &traitPosition,
	     lDT[0], 275, lDT[1], 275, lDT[2], 275, lDT[3], 275, lDT[4], 275,
	     lDT[5], 275);
  tpl_pack (tn, 0);
  tpl_dump (tn, TPL_FILE, fileName);
  tpl_free (tn);

  return 0;
}

int
restoreAlternative (char *pedigree, int chromosome, double traitPosition,
		    double **lDT)
{
  tpl_node *tn;
  FILE *file;
  char *checkPedigree;
  int checkChromosome;
  double checkTraitPosition;

  //  fprintf(stderr, "in restoreAlternative for pedigree %s, chromosome %d, trait position %d\n",
  //      pedigree, chromosome, (int) traitPosition);
  sprintf (fileName, alternativeFileFormat, resultsprefix, pedigree,
	   chromosome, traitPosition);
  tn =
    tpl_map (alternativeTPLFormat, &checkPedigree, &checkChromosome,
	     &checkTraitPosition, lDT[0], 275, lDT[1], 275, lDT[2], 275,
	     lDT[3], 275, lDT[4], 275, lDT[5], 275);
  if ((file = fopen (fileName, "r"))) {
    fclose (file);
    tpl_load (tn, TPL_FILE, fileName);
    tpl_unpack (tn, 0);
    if ((strcmp (pedigree, checkPedigree) == 0)
	&& (chromosome == checkChromosome)
	&& (traitPosition == checkTraitPosition)) {
      tpl_free (tn);
      //      fprintf(stderr, "restoreAlternative successfully restored alternative data\n");
      return 1;
    } else {
      fprintf (stderr,
	       "restoreAlternative check of pedigree/chromosome/traitPosition %s/%d/%G vs %s/%d/%G failed, exiting!\n",
	       pedigree, chromosome, traitPosition, checkPedigree,
	       checkChromosome, checkTraitPosition);
      exit (1);
    }
  }
  //  fprintf(stderr, "restoreAlternative found no alternative data\n");
  tpl_free (tn);
  return 0;
}


#ifdef TEST_MAIN

#include <math.h>

int
main ()
{
  int i, j;
  double **lDT;

  if ((lDT = (double **) malloc (6 * sizeof (double *))) == NULL) {
    fprintf (stderr, "malloc of 1st dimension failed!\n");
    exit (1);
  }
  for (i = 0; i < 6; i++) {
    if ((lDT[i] = (double *) malloc (275 * sizeof (double))) == NULL) {
      fprintf (stderr, "malloc of 2nd dimension failed\n");
      exit (1);
    }
    for (j = 0; j < 275; j++)
      lDT[i][j] = (double) i + sqrt (j + 1);
  }
  saveTrait ("test", lDT);
  free (lDT);
  lDT = restoreTrait ("test");
  dump_lDT (lDT);
  return 0;
}

#endif
