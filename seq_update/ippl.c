#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "ippl.h"

st_pplchr *pplchrs = NULL;
int numpplchrs = 0;
char *pplinfile;

#define IPPL_CHR_COL    1
#define IPPL_POS_COL    2
#define IPPL_PPL_COL    3

#define BUFFLEN 256

void read_ppls (char *filename)
{
  int numcols=0, *datacols=NULL, lineno=0, foundcols=0, colno;
  double pos, ppl;
  char buff[BUFFLEN], chr[8], *token, *pb, *endptr;
  FILE *fp;
  
  pplinfile = filename;
  if ((fp = fopen (pplinfile, "r")) == NULL) {
    fprintf (stderr, "open '%s' failed, %s\n", pplinfile, strerror (errno));
    exit (-1);
  }

  /* Not requiring that the PPL file have a version line, since the multipoint PPLs might
   * have been calculated by something other than Kelvin.
   */
  while (1) {
    if (fgets (buff, BUFFLEN, fp) == NULL) {
      if (feof (fp))
	fprintf (stderr, "pplfile '%s' ends unexpectedly at line %d\n", pplinfile, lineno);
      else 
	fprintf (stderr, "error reading '%s' at line %d, %s\n", pplinfile, lineno,
		 strerror (errno));
      exit (-1);
    }
    lineno++;
    if (strcasecmp (buff, "# Version") != 0)
      break;
  }

  /* We'l try to be flexible about file format, so long as there are column headers that
   * allow us to identify the columns. 
   */

  token = strtok_r (buff, " \t\n", &pb);
  while ((token != NULL) &&
	 (foundcols != ((1 << IPPL_CHR_COL) | (1 << IPPL_POS_COL) | (1 << IPPL_PPL_COL)))) {
    numcols++;
    if ((datacols = (int *) realloc (datacols, sizeof (int) * numcols)) == NULL) {
      fprintf (stderr, "realloc datacols failed, %s\n", strerror (errno));
      exit (-1);
    }
    
    if ((strcasecmp (token, "Chr") == 0) || (strcasecmp (token, "Chromosome") == 0)) {
      /* The chromosome column */
      datacols[numcols - 1] = IPPL_CHR_COL;
      foundcols |= 1 << IPPL_CHR_COL;

    } else if ((strcasecmp (token, "Pos") == 0) || (strcasecmp (token, "Position") == 0)) {
      /* The position column */
      datacols[numcols - 1] = IPPL_POS_COL;
      foundcols |= 1 << IPPL_POS_COL;
      
    } else if (strcasecmp (token, "PPL") == 0) {
      datacols[numcols - 1] = IPPL_PPL_COL;
      foundcols |= 1 << IPPL_PPL_COL;
      
    } else {
      /* Something else */
      datacols[numcols - 1] = 0;
    }
    
    token = strtok_r (NULL, " \t\n", &pb);
  }
  if ((foundcols & (1 << IPPL_CHR_COL)) == 0) {
    fprintf (stderr, "pplfile '%s' contains no chromosome column\n", pplinfile);
    exit (-1);
  } else if ((foundcols & (1 << IPPL_POS_COL)) == 0) {
    fprintf (stderr, "pplfile '%s' contains no position column\n", pplinfile);
    exit (-1);
  } else if ((foundcols & (1 << IPPL_PPL_COL)) == 0) {
    fprintf (stderr, "pplfile '%s' contains no PPL column\n", pplinfile);
    exit (-1);
  }

  while (fgets (buff, BUFFLEN, fp) != NULL) {
    lineno++;
    token = strtok_r (buff, " \t\n", &pb);
    for (colno = 0; colno < numcols; colno++) {
      if (token == NULL) {
	fprintf (stderr, "short line in '%s', line %d\n", pplinfile, lineno);
	exit (-1);
      }
      switch (datacols[colno]) {
      case IPPL_CHR_COL:
	/* Sure, I could handle this better, but seriously, how many characters do
	 * you need to specify the chromosome? 
	 */
	if (strlen (token) > 7) {
	  fprintf (stderr, "chromosome identifier '%s' in pplfile '%s' is too long\n", 
		   token, pplinfile);
	    exit (-1);
	}
	endptr = strcpy (chr, token);
	break;
      case IPPL_POS_COL:
	pos = strtod (token, &endptr);
	break;
      case IPPL_PPL_COL:
	ppl = strtod (token, &endptr);
	break;
      }
      if (token == endptr) {
	fprintf (stderr, "illegal data in line %d in '%s'\n", lineno, pplinfile);
	exit (-1);
      }
      token = strtok_r (NULL, " \t\n", &pb);
    }
    insert_ppls (chr, pos, ppl);
  }
  if (! feof (fp)) {
    fprintf (stderr, "error reading '%s' at line %d, %s\n", pplinfile, lineno, strerror (errno));
    exit (-1);
  }
  fclose (fp);
  free (datacols);
  return;
}


void insert_ppls (char *chr, double pos, double ppl)
{
  int va;
  st_pplchr *ptr = NULL;

  for (va = 0; va < numpplchrs; va++) {
    if (strcmp (pplchrs[va].chr, chr) == 0) {
      ptr = &pplchrs[va];
      break;
    }
  }
  if (ptr == NULL) {
    if ((pplchrs = (st_pplchr *) realloc (pplchrs, sizeof (st_pplchr) * (numpplchrs + 1))) == NULL) {
      fprintf (stderr, "realloc pplchrs failed, %s\n", strerror (errno));
      exit (-1);
    }
    ptr = &pplchrs[numpplchrs];
    numpplchrs++;
    memset (ptr, 0, sizeof (st_pplchr));
    strcpy (ptr->chr, chr);
  }

  if ((ptr->numpos != 0) && (ptr->pos[ptr->numpos-1] > pos)) {
    fprintf (stderr, "position '%f' out of order in '%s'\n", pos, pplinfile);
    exit (-1);
  }
  if (((ptr->pos = (double *) realloc (ptr->pos, sizeof (double) * ptr->numpos+1)) == NULL) || 
      ((ptr->ppl = (double *) realloc (ptr->ppl, sizeof (double) * ptr->numpos+1)) == NULL)) {
    fprintf (stderr, "realloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  ptr->pos[ptr->numpos] = pos;
  ptr->ppl[ptr->numpos] = ppl;
  ptr->numpos++;
  return;
}


double get_ippl (char *chr, double pos)
{
  int va;
  double ippl;
  st_pplchr *ptr = NULL;

  for (va = 0; va < numpplchrs; va++) {
    if (strcmp (pplchrs[va].chr, chr) == 0) {
      ptr = &pplchrs[va];
      break;
    }
  }
  if (ptr == NULL) {
    fprintf (stderr, "'%s' doesn't contain data for chromosome '%s'\n", pplinfile, chr);
    exit (-1);
  }
  if (pos < ptr->pos[0])
    return (ptr->ppl[0]);
  if (pos > ptr->pos[ptr->numpos-1])
    return (ptr->ppl[ptr->numpos-1]);
  for (va = 1; va < ptr->numpos; va++) {
    if ((ptr->pos[va-1] <= pos) && (pos <= ptr->pos[va])) {
      
      ippl = ptr->ppl[va-1] + (ptr->ppl[va] - ptr->ppl[va-1]) *
	((pos - ptr->pos[va-1]) / (ptr->pos[va] - ptr->pos[va-1]));
      return (ippl);
    }
  }
}


void free_ppls ()
{
  int va;

  for (va = 0; va < numpplchrs; va++) {
    if (pplchrs[va].numpos > 0) {
      free (pplchrs[va].pos);
      free (pplchrs[va].ppl);
    }
  }
  if (numpplchrs > 0)
    free (pplchrs);

  return;
}
