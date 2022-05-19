/**
@file saveResults.c

  Routines to save and restore trait, marker and alternative likelhood
  result grids from a 6 disease gene frequency by 275 penetrance
  analysis. Eventually this will be a variable-sized grid.  The data
  being passed is both self-referential and not contiguous in memory,
  so it cannot be serialized all in one shot.

  Uses TPL serialization, which is open source under the Berkeley license.

  Author: Bill Valentine-Cooper

  Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
  This program is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.
  You should have received a copy of the GNU General Public License along
  with this program. If not, see <https://www.gnu.org/licenses/>.

*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "kelvin.h"
#include "saveResults.h"
#include "utils/tpl.h"
#include "config/config.h"

char pathName[256];     ///< Where the path to saved results will be built. Should malloc.
char fileName[256];     ///< Filename for saved results files. Should malloc.

/**

  Dump a 6x275 grid of doubles. Used to prove that we are actually
  saving what we're passed, and restoring what we saved.


*/
void dump_lDT (double **lDT     ///< Pointer to list of 6x275 pointers to doubles.
  )
{
  int i, j;

  for (i = 0; i < 6; i++)
    for (j = 0; j < 275; j++)
      fprintf (stderr, "%G ", lDT[i][j]);
  fprintf (stderr, "\n");
  return;
}

char *traitTPLFormat = "sf#f#f#f#f#f#"; ///< String and six fixed vectors of doubles for trait
char *traitFileFormat = "%sped-%s_trait.tpl";   ///< Trait filename components
/**

  Save a 6x275 grid of doubles of trait information along with descriptive attributes.
  Use either:

  <pre>
  <results-prefix>trait-23/ped-<pid>/ped-<pid>_trait.tpl
  for chromosome 23, or
  <results-prefix>trait/ped-<pid>/ped-<pid>_trait.tpl
  for other chromosomes

  where:

  <results-prefix> is the parameter to the SR directive in the configuration
      file, or "./" by default.
  <pid> is the pedigree ID string
  </pre>

  Always returns a 0.

*/
int saveTrait (int chr23Flag,   ///< Flag is 1 if doing chromosome 23, 0 otherwise.
               char *pedigree,  ///< Name of the pedigree to use in path and filename.
               double **lDT     ///< Pointer to list of 6x275 pointers to doubles of likelihoods.
  )
{
  tpl_node *tn;

  mkdir (modelOptions->resultsprefix, S_IRWXU | S_IRWXG | S_IROTH);
  if (chr23Flag) {
    sprintf (pathName, "%strait-23/", modelOptions->resultsprefix);
    mkdir (pathName, S_IRWXU | S_IRWXG | S_IROTH);
    sprintf (pathName, "%strait-23/ped-%s/", modelOptions->resultsprefix, pedigree);
  } else {
    sprintf (pathName, "%strait/", modelOptions->resultsprefix);
    mkdir (pathName, S_IRWXU | S_IRWXG | S_IROTH);
    sprintf (pathName, "%strait/ped-%s/", modelOptions->resultsprefix, pedigree);
  }
  mkdir (pathName, S_IRWXU | S_IRWXG | S_IROTH);
  sprintf (fileName, traitFileFormat, pathName, pedigree);
  tn = tpl_map (traitTPLFormat, &pedigree, lDT[0], 275, lDT[1], 275, lDT[2], 275, lDT[3], 275, lDT[4], 275, lDT[5], 275);
  tpl_pack (tn, 0);
  tpl_dump (tn, TPL_FILE, fileName);
  tpl_free (tn);
  return 0;
}

/**

  Restore any saved trait results matching the parameters
  provided. Parameters are identical to saveTrait.

  Return a 0 if no matching saved trait results are found. Internal
  consistency is checked as well, i.e. if the parameter values do not
  match file contents, the program will exit.

*/
int restoreTrait (int chr23Flag, char *pedigree, double **lDT)
{
  tpl_node *tn;
  FILE *file;
  char *checkPedigree;

  if (chr23Flag)
    sprintf (pathName, "%strait-23/ped-%s/", modelOptions->resultsprefix, pedigree);
  else
    sprintf (pathName, "%strait/ped-%s/", modelOptions->resultsprefix, pedigree);
  sprintf (fileName, traitFileFormat, pathName, pedigree);
  tn = tpl_map (traitTPLFormat, &checkPedigree, lDT[0], 275, lDT[1], 275, lDT[2], 275, lDT[3], 275, lDT[4], 275, lDT[5], 275);
  if ((file = fopen (fileName, "r"))) {
    fclose (file);
    tpl_load (tn, TPL_FILE, fileName);
    tpl_unpack (tn, 0);
    if (strcmp (pedigree, checkPedigree) == 0) {
      tpl_free (tn);
      return 1;
    } else {
      fprintf (stderr,
               "restoreTrait check of requested pedigree of %s vs restored of %s failed, exiting!\n", pedigree, checkPedigree);
      exit (1);
    }
  }
  tpl_free (tn);
  return 0;
}

char *markerTPLFormat = "siiA(s)f";     ///< String, two ints, array of string and a single float
char *markerFileFormat = "%schr-%d_ped-%s_";    ///< Marker filename components
/**

  Save a 6x275 grid of doubles of marker information along with descriptive attributes.
  Use path and name:

  <pre>
  <results-prefix>chr-<chr>/ped-<pid>/chr-<chr>_ped-<pid>_<mk1>_<mk2>{_<mkn>...}_marker.tpl

  where:

  <results-prefix> is the parameter to the SR directive in the configuration file, or "./" by default.
  <chr> is the chromosome number
  <pid> is the pedigree ID string
  <mk1> is the name of the first marker
  <mk2> is the name of the second marker
  <mkn> is the name of the nth marker
  <pre>

  Always returns a 0.

*/
int saveMarker (char *pedigree, ///< Name of pedigree to use in path and filename.
                int chromosome, ///< Chromosome number to use in path and filename.
                int markerCount,        ///< Number of markers in analysis, currently just 2.
                char **markerNames,     ///< List of pointers to marker names.
                double *mDT     ///< Pointer to double of marker likelihood.
  )
{
  tpl_node *tn;
  int i;
  char *markerName;

  mkdir (modelOptions->resultsprefix, S_IRWXU | S_IRWXG | S_IROTH);
  sprintf (pathName, "%schr-%d/", modelOptions->resultsprefix, chromosome);
  mkdir (pathName, S_IRWXU | S_IRWXG | S_IROTH);
  sprintf (pathName, "%schr-%d/ped-%s/", modelOptions->resultsprefix, chromosome, pedigree);
  mkdir (pathName, S_IRWXU | S_IRWXG | S_IROTH);
  sprintf (fileName, markerFileFormat, pathName, chromosome, pedigree);
  for (i = 0; i < markerCount; i++) {
    strcat (fileName, markerNames[i]);
    strcat (fileName, "_");
  }
  strcat (fileName, "marker.tpl");
  tn = tpl_map (markerTPLFormat, &pedigree, &chromosome, &markerCount, &markerName, mDT);
  tpl_pack (tn, 0);

  for (i = 0; i < markerCount; i++) {
    markerName = markerNames[i];
    tpl_pack (tn, 1);
  }

  tpl_dump (tn, TPL_FILE, fileName);
  tpl_free (tn);

  return 0;
}


/**

  Restore any saved marker results matching the parameters
  provided. Parameters are identical to saveTrait.

  Return a 0 if no matching saved marker results are found. Internal
  consistency is checked as well, i.e. if the parameter values do not
  match file contents, the program will exit.

*/
int restoreMarker (char *pedigree, int chromosome, int markerCount, char **markerNames, double *mDT)
{
  int i, checkChromosome, checkMarkerCount;
  tpl_node *tn;
  FILE *file;
  char *checkPedigree;
  char *markerName;

  sprintf (pathName, "%schr-%d/ped-%s/", modelOptions->resultsprefix, chromosome, pedigree);
  sprintf (fileName, markerFileFormat, pathName, chromosome, pedigree);
  for (i = 0; i < markerCount; i++) {
    strcat (fileName, markerNames[i]);
    strcat (fileName, "_");
  }
  strcat (fileName, "marker.tpl");
  tn = tpl_map (markerTPLFormat, &checkPedigree, &checkChromosome, &checkMarkerCount, &markerName, mDT);
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
      return 1;
    } else {
      fprintf (stderr,
               "restoreMarker check of pedigree/marker/count of %s/%d/%d vs %s/%d/%d failed, exiting!\n",
               pedigree, chromosome, markerCount, checkPedigree, checkChromosome, checkMarkerCount);
      exit (1);
    }
  }
  tpl_free (tn);
  return 0;
}

char *alternativeTPLFormat = "siff#f#f#f#f#f#"; /* String and six fixed vectors of doubles */
char *alternativeFileFormat = "%schr-%d_ped-%s_pos-%G_alternative.tpl";
/**

  Save a 6x275 grid of doubles of alternative likelihood information along with descriptive
  attributes. Use:

  <pre>
  <results-prefix>chr-<chr>/ped-<pid>/chr-<chr>_ped-<pid>/chr-<chr>_ped-<pid>_pos-<trt>_alternative.tpl

  where:

  <results-prefix> is the parameter to the SR directive in the configuration file, or "./" by default.
  <chr> is the chromosome number
  <pid> is the pedigree ID string
  <trt> is the trait position number
  </pre>

  Always returns a 0.

*/
int saveAlternative (char *pedigree,    ///< Name of pedigree to use in path and filename.
                     int chromosome,    ///< Chromosome number to use in path and filename.
                     double traitPosition,      ///< Trait position to use in path and filename.
                     double **lDT       ///< Pointer to list of pointers to doubles of alternative likelihood.
  )
{
  tpl_node *tn;

  mkdir (modelOptions->resultsprefix, S_IRWXU | S_IRWXG | S_IROTH);
  sprintf (pathName, "%schr-%d/", modelOptions->resultsprefix, chromosome);
  mkdir (pathName, S_IRWXU | S_IRWXG | S_IROTH);
  sprintf (pathName, "%schr-%d/ped-%s/", modelOptions->resultsprefix, chromosome, pedigree);
  mkdir (pathName, S_IRWXU | S_IRWXG | S_IROTH);
  sprintf (fileName, alternativeFileFormat, pathName, chromosome, pedigree, traitPosition);
  tn = tpl_map (alternativeTPLFormat, &pedigree, &chromosome, &traitPosition,
                lDT[0], 275, lDT[1], 275, lDT[2], 275, lDT[3], 275, lDT[4], 275, lDT[5], 275);
  tpl_pack (tn, 0);
  tpl_dump (tn, TPL_FILE, fileName);
  tpl_free (tn);

  return 0;
}

/**

  Restore any saved alternative results matching the parameters
  provided. Parameters are identical to saveTrait.

  Return a 0 if no matching saved alternative results are found. Internal
  consistency is checked as well, i.e. if the parameter values do not
  match file contents, the program will exit.

*/
int restoreAlternative (char *pedigree, int chromosome, double traitPosition, double **lDT)
{
  tpl_node *tn;
  FILE *file;
  char *checkPedigree;
  int checkChromosome;
  double checkTraitPosition;

  sprintf (pathName, "%schr-%d/ped-%s/", modelOptions->resultsprefix, chromosome, pedigree);
  sprintf (fileName, alternativeFileFormat, pathName, chromosome, pedigree, traitPosition);
  tn = tpl_map (alternativeTPLFormat, &checkPedigree, &checkChromosome,
                &checkTraitPosition, lDT[0], 275, lDT[1], 275, lDT[2], 275, lDT[3], 275, lDT[4], 275, lDT[5], 275);
  if ((file = fopen (fileName, "r"))) {
    fclose (file);
    tpl_load (tn, TPL_FILE, fileName);
    tpl_unpack (tn, 0);
    if ((strcmp (pedigree, checkPedigree) == 0)
        && (chromosome == checkChromosome)
        && (traitPosition == checkTraitPosition)) {
      tpl_free (tn);
      return 1;
    } else {
      fprintf (stderr,
               "restoreAlternative check of pedigree/chromosome/traitPosition %s/%d/%G vs %s/%d/%G failed, exiting!\n",
               pedigree, chromosome, traitPosition, checkPedigree, checkChromosome, checkTraitPosition);
      exit (1);
    }
  }
  tpl_free (tn);
  return 0;
}
