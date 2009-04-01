/*
 * calc_updated_ppl - Sequentially update Kelvin-format BR files, and
 *                    calculate PPL and various LD related PPL-like statistics
 *
 * John Burian - john.burian@nationwidechildrens.org
 *
 * Copyright 2008, The Research Institute at Nationwide Children's Hospital
 * All rights reserved. Permission is granted to use this software for
 * non-profit educational purposes only.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <getopt.h>

/* Kelvin header files */
#include <kelvin.h>
#include <integrationLocals.h>
#include <multidim.h>
#include <ippl.h>
#include <map.h>

/* TO DO:
 * - sex-specific theta cutoffs
 * - mechanism for detecting missing D'/Theta combinations
 */

#define DPRIME_COL 1
#define THETA_COL  2
#define LR_COL     3
#define POS_COL    4
#define CHR_COL    5
#define PPL_COL    6

#define METH_OLD   1
#define METH_NEW   2

#define BUFFLEN 256
/* Default prior for LD statistic calculations */
#define DEFAULT_LDPRIOR 0.02
/* the multipoint PPL for a BR of 0.214 */
#define MIN_PRIOR 7.8528124097619317e-03

typedef struct {
  char name1[MAX_MAP_NAME_LEN],
    name2[MAX_MAP_NAME_LEN],
    chr[MAX_MAP_CHR_LEN];
  double pos;
  int num;
} st_brmarker;

typedef struct {
  char *name;
  FILE *fp;
  int lineno,
    numcols,
    numdprimes,
    numthetas,
    no_ld,
    holey_grid,
    two_point,
    eof,
    datacolsize,
    *datacols;
  st_brmarker curmarker;
} st_brfile;

typedef struct {
  int dprimesize,
    thetasize;
  double lr,
    *dprimes,
    *thetas;
} st_data;

typedef LDVals st_ldvals;
/* LDVals looks like this :
 *
typedef struct {
  double ld_small_theta,
    ld_big_theta,
    ld_unlinked,
    le_small_theta,
    le_big_theta,
    le_unlinked;
} st_ldvals;
*/


/* Defaults for varaiables that define the nature of the data */
int sexspecific = 0;    /* -s or --sexspecific: sex-specific thetas */
int multipoint = 0;     /* -m or --multipoint */
int okelvin = 0;        /* -o or --okelvin : expect original (fixed-grid) kelvin BR data */
int relax = 0;          /* -r or --relax: don't compare marker names between files */
int allstats = 0;       /* -a or --allstats: print all LD and iPPL stats */
int quiet = 0;          /* -q or --quiet : supress non-fatal warnings */

/* Defaults for variables used in calculation of the statistics */
double prior = 0.02;    /* -p or --prior: prior probability */
double weight = 0.95;   /* -w or --weight: weighting factor */
double cutoff = 0.05;   /* -c or --cutoff: theta cutoff */

/* For input BRs with fake Chr/cM info, force the output to use Chr/cM from the optional map */
int forcemap = 0;

/* For optional output files, the filehandles are global */
FILE *partout = NULL;   /* --partout: file to which to write partial results */
FILE *bfout = NULL;     /* --bfout: file to which bayes factors will be written */
FILE *sixout = NULL;    /* --sixout: file to which six-region (ldval) values will be written */

/* For option input files, the filenames are global */
char *pplinfile=NULL;   /* --pplin: file from which to read PPLs for use as position-
	  		            specific priors */
char *mapinfile=NULL;   /* --mapin: mapfile, to allow updating across files containing sets of
 			            markers that only partially overlap */


/* Globally available for error messages */
char *current = "0.36.1";
char *pname;
int verbose = 0;

void kelvin_twopoint (st_brfile *brfiles, int numbrfiles);
void kelvin_multipoint (st_brfile *brfiles, int numbrfiles);
void dkelvin_twopoint (st_brfile *brfiles, int numbrfiles);
void print_twopoint_headers (int no_ld);
void print_twopoint_stats (int no_ld, st_brmarker *marker, st_ldvals *ldval);
void do_first_pass (st_brfile *brfile, st_multidim *dprimes, st_multidim *thetas, st_data *data);
void do_dkelvin_first_pass (st_brfile *brfile, st_data *data);
/*
double calc_ppl_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr);
double calc_ppl_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr);
*/
void calc_ldvals_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr, st_ldvals *ldval);
void calc_ldvals_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr, st_ldvals *ldval);
double calc_upd_ppl (st_ldvals *ldval);
double calc_upd_ppl_allowing_ld (st_ldvals *ldval, double ldprior);
double calc_upd_ppld_given_linkage (st_ldvals *ldval, double ldprior);
double calc_upd_ppld_allowing_l (st_ldvals *ldval, double ldprior);

int parse_command_line (int argc, char **argv);
void usage ();
double validate_double_arg (char *arg, char *optname);
void open_brfile (st_brfile *brfile);
void get_next_marker (st_brfile *brfile, st_data *data);
int get_marker_line (st_brfile *brfile);
int get_header_line (st_brfile *brfile);
int get_data_line (st_brfile *brfile, st_data *data);
void print_partial_header (st_brfile *brfile);
void print_partial_data (st_brfile *brfile, st_data *data);
int compare_markers (st_brmarker *m1, st_brmarker *m2);
void compare_headers (st_brfile *f1, st_brfile *f2);
int compare_positions (st_brmarker *m1, st_brmarker *m2);
void compare_samples (st_data *d, int sampleno, st_brfile *brfile);


int main (int argc, char **argv)
{
  int argidx, va, numbrfiles;
  st_brfile *brfiles;

  pname = argv[0];
  argidx = parse_command_line (argc, argv);

  if (pplinfile != NULL)
    read_ppls (pplinfile);
  if (mapinfile != NULL)
    read_map (mapinfile);

  if (argidx >= argc) {
    fprintf (stderr, "missing file name\n");
    exit (-1);
  }
  numbrfiles = argc - argidx;
  if ((brfiles = calloc (numbrfiles, sizeof (st_brfile))) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < numbrfiles; va++) {
    brfiles[va].name = argv[argidx + va];
    open_brfile (&brfiles[va]);
  }
  printf ("# Version V%s\n", current);
  if (partout != NULL)
    fprintf (partout, "# Version V%s\n", current);
  
  if (multipoint) {
    kelvin_multipoint (brfiles, numbrfiles);
  } else if (okelvin) {
    kelvin_twopoint (brfiles, numbrfiles);
  } else {
    dkelvin_twopoint (brfiles, numbrfiles);
  }
  
  for (va = 0; va < numbrfiles; va++) {
    fclose (brfiles[va].fp);
  }
  free (brfiles);
  if (pplinfile != NULL)
    free_ppls ();
  if (mapinfile != NULL)
    free_map ();
  if (partout != NULL)
    fclose (partout);
  if (bfout != NULL)
    fclose (bfout);
  exit (0);
}


void kelvin_twopoint (st_brfile *brfiles, int numbrfiles)
{
  int fileno, alldone=0, didx, thidx, numcurrent, dsize=0, thsize=0, ret;
  double **lr=NULL, ldprior, ldstat;
  st_brmarker next_marker;
  st_mapmarker *mapptr;
  st_data data;
  st_multidim dprimes, thetas;
  st_brfile **current;
  st_ldvals ldval;

  memset (&next_marker, 0, sizeof (st_brmarker));
  memset (&data, 0, sizeof (st_data));
  memset (&dprimes, 0, sizeof (st_multidim));
  memset (&thetas, 0, sizeof (st_multidim));
  memset (&ldval, 0, sizeof (st_ldvals));
  if ((current = malloc (sizeof (st_brfile *) * numbrfiles)) == NULL) {
    fprintf (stderr, "malloc current failed, %s\n", strerror (errno));
    exit (-1);
  }

  for (fileno = 0; fileno < numbrfiles; fileno++) {
    get_next_marker (&brfiles[fileno], &data);
    if (verbose >= 2)
      printf ("first marker in %s is %s\n", brfiles[fileno].name, brfiles[fileno].curmarker.name2);
  }
  print_twopoint_headers (brfiles[0].no_ld);

  while (1) {
    if (mapinfile != NULL) {
      if ((mapptr = next_mapmarker (NULL)) != NULL) {
	strcpy (next_marker.chr, mapptr->chr);
	strcpy (next_marker.name2, mapptr->name);
      } else {
	break;
      }
    } else {
      memcpy (&next_marker, &brfiles[0].curmarker, sizeof (st_brmarker));
    }
    
    if (verbose >= 2)
      printf ("first marker to update is %s\n", next_marker.name2);
    
    numcurrent = 0;
    alldone = 1;
    for (fileno = 0; fileno < numbrfiles; fileno++) {
      if (brfiles[fileno].eof) {
	if (verbose >= 3)
	  printf ("%s is at EOF\n", brfiles[fileno].name);
	continue;
      }
      alldone = 0;
      if (compare_markers (&next_marker, &brfiles[fileno].curmarker) == -1) {
	if (mapinfile != NULL)
	  continue;
	fprintf (stderr, "marker mismatch at line %d in '%s', expecting %s, found %s\n",
		 brfiles[fileno].lineno, brfiles[fileno].name, next_marker.name2,
		 brfiles[fileno].curmarker.name2);
	exit (-1);
      }
      current[numcurrent++] = &brfiles[fileno];
      compare_headers (current[0], &brfiles[fileno]);
    }
    if (alldone)
      break;
    if (numcurrent == 0) {
      if (! quiet)
	fprintf (stderr, "warning: marker %s from %s doesn't appear in any BR files, skipping\n",
		 next_marker.name2, mapinfile);
      continue;
    }
    if (verbose >= 1) {
      printf ("current marker is '%s', using BR file%s %s", next_marker.name2,
	      (numcurrent > 1) ? "s" : "", current[0]->name);
      for (fileno = 1; fileno < numcurrent; fileno++) 
	printf (", %s", current[fileno]->name);
      printf ("\n");
    }
    
    if (partout != NULL)
      print_partial_header (current[0]);
    
    do_first_pass (current[0], &dprimes, &thetas, &data);
    if (dsize < dprimes.totalelems) {
      if ((lr = realloc (lr, sizeof (double *) * dprimes.totalelems)) == NULL) {
	fprintf (stderr, "malloc failed, %s\n", strerror (errno));
	exit (-1);	
      }
      for ( ; dsize < dprimes.totalelems; dsize++)
	lr[dsize] = NULL;
    }
    for (didx = 0; didx < dsize; didx++) {
      if ((thsize < thetas.totalelems) || (lr[didx] == NULL)) {
	if ((lr[didx] = realloc (lr[didx], sizeof (double) * thetas.totalelems)) == NULL) {
	  fprintf (stderr, "malloc failed, %s\n", strerror (errno));
	  exit (-1);
	}
      }
    }
    if (thsize < thetas.totalelems)
      thsize = thetas.totalelems;
    
    for (didx = 0; didx < dprimes.totalelems; didx++) {
      for (thidx = 0; thidx < thetas.totalelems; thidx++) {
	lr[didx][thidx] = DBL_MAX;
      }
    }
    
    while (get_data_line (current[0], &data) == 1) {
      if ((didx = multi_find (&dprimes, data.dprimes, current[0]->numdprimes)) == -1) {
	fprintf (stderr, "unexpected dprime(s) in '%s' at line %d\n", current[0]->name,
		 current[0]->lineno);
	exit (-1);
      }
      if ((thidx = multi_find (&thetas, data.thetas, current[0]->numthetas)) == -1) {
	fprintf (stderr, "unexpected thetas in '%s' at line %d\n", current[0]->name,
		 current[0]->lineno);
	exit (-1);
      }
      if (lr[didx][thidx] != DBL_MAX) {
	fprintf (stderr, "duplicate dprime/theta combination at line %d in '%s'\n",
		 current[0]->lineno, current[0]->name);
	exit (-1);
      }
      lr[didx][thidx] = data.lr;
      
      for (fileno = 1; fileno < numcurrent; fileno++) {
	if ((ret = get_data_line (current[fileno], &data)) != 1) {
	  fprintf (stderr, "expected data in '%s' at line %d, found %s\n",
		   current[fileno]->name, current[fileno]->lineno,
		   (ret == 0) ? "end-of-file" : "marker line");
	  exit (-1);
	}
	if ((ret = multi_find (&dprimes, data.dprimes, current[fileno]->numdprimes)) != didx) {
	  fprintf (stderr, "%s dprime(s) in '%s' at line %d\n", (ret == -1) ? "unexpected" :
		   "misordered", current[fileno]->name, current[fileno]->lineno);
	  exit (-1);
	}
	if ((ret = multi_find (&thetas, data.thetas, current[fileno]->numthetas)) != thidx) {
	  fprintf (stderr, "unexpected thetas in '%s' at line %d\n", (ret == -1) ? "unexpected" :
		   "misordered", current[fileno]->name, current[fileno]->lineno);
	  exit (-1);
	}
	if (lr[didx][thidx] == DBL_MAX) {
	  fprintf (stderr, "dprime/theta combination at line %d in '%s' does not appear in '%s' \n",
		   current[fileno]->lineno, current[fileno]->name, current[0]->name);
	  exit (-1);
	}
	lr[didx][thidx] *= data.lr;
      }
      if (partout != NULL) {
	data.lr = lr[didx][thidx];
	print_partial_data (current[0], &data);
      }
    }

    if (! sexspecific)
      calc_ldvals_sexavg (&dprimes, &thetas, lr, &ldval);
    else 
      calc_ldvals_sexspc (&dprimes, &thetas, lr, &ldval);

    if (forcemap) {
      strcpy (current[0]->curmarker.chr, mapptr->chr);
      current[0]->curmarker.pos = mapptr->avgpos;
    }
    print_twopoint_stats (current[0]->no_ld, &(current[0]->curmarker), &ldval);
    for (fileno = 0; fileno < numcurrent; fileno++) {
      get_next_marker (current[fileno], &data);
    }
  }

  for (didx = 0; didx < dsize; didx++) {
    if (lr[didx] != NULL)
      free (lr[didx]);
  }
  if (lr != NULL)
    free (lr);
  multi_free (&dprimes);
  multi_free (&thetas);
  free (current);
  if (data.dprimesize > 0)
    free (data.dprimes);
  if (data.thetasize > 0)
    free (data.thetas);
  for (fileno = 0; fileno < numbrfiles; fileno++) {
    if (brfiles[fileno].datacolsize > 0)
      free (brfiles[fileno].datacols);
  }
  return;
}


void kelvin_multipoint (st_brfile *brfiles, int numbrfiles)
{
  int fileno, ret, alldone, curno;
  st_brmarker *marker;
  st_data data;
  double lr, *lrs, ppl;

  memset (&data, 0, sizeof (st_data));
  if ((lrs = malloc (sizeof (double) * numbrfiles)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  
  get_header_line (&brfiles[0]);
  if ((ret = get_data_line (&brfiles[0], &data)) != 1) {
    fprintf (stderr, "expected data in '%s' at line %d, found %s\n", brfiles[0].name,
	     brfiles[0].lineno, (ret == 0) ? "end-of-file" : "marker line");
    exit (-1);
  }
  lrs[0] = data.lr;
  for (fileno = 1; fileno < numbrfiles; fileno++) {
    get_header_line (&brfiles[fileno]);
    compare_headers (&brfiles[0], &brfiles[fileno]);
    if ((ret = get_data_line (&brfiles[fileno], &data)) != 1) {
      fprintf (stderr, "expected data in '%s' at line %d, found %s\n", brfiles[fileno].name,
	       brfiles[fileno].lineno, (ret == 0) ? "end-of-file" : "marker line");
      exit (-1);
    }
    lrs[fileno] = data.lr;
  }
  printf ("Chr Position PPL BayesRatio\n");

  while (1) {
    alldone = 1;
    curno = -1;
    lr = 1;
    for (fileno = 0; fileno < numbrfiles; fileno++) {
      if (brfiles[fileno].eof)
	continue;
      alldone = 0;
      if (curno == -1)
	curno = fileno;
      else if (compare_positions (&brfiles[fileno].curmarker, &brfiles[curno].curmarker) < 0)
	curno = fileno;
    }
    if (alldone)
      break;
    for (fileno = 0; fileno < numbrfiles; fileno++) {
      if (brfiles[fileno].eof)
	continue;
      if (compare_positions (&brfiles[fileno].curmarker, &brfiles[curno].curmarker) != 0)
	continue;
      lr *= lrs[fileno];
      if (fileno != curno) {
	if (get_data_line (&brfiles[fileno], &data) == 1)
	  lrs[fileno] = data.lr;
      }
    }

    marker = &brfiles[curno].curmarker;
    if ((lr < 0.214) || ((ppl = (lr * lr) / (-5.77 + (54 * lr) + (lr * lr))) < 0.0))
      ppl = 0.0;
    printf ("%s %.4f %.3f %.6e\n", marker->chr, marker->pos, ppl, lr);
    
    if (get_data_line (&brfiles[curno], &data) == 1)
      lrs[curno] = data.lr;
  }
  free (lrs);
  return;
}


void dkelvin_twopoint (st_brfile *brfiles, int numbrfiles)
{
  int fileno, numcurrent, sampleno, alldone, ret;
  double ldprior, ldstat;
  st_brmarker next_marker, *marker;
  st_mapmarker *mapptr;
  st_data data_0, data_n;
  st_ldvals ldval;
  st_brfile **current;

  memset (&next_marker, 0, sizeof (st_brmarker));
  memset (&data_0, 0, sizeof (st_data));
  memset (&data_n, 0, sizeof (st_data));
  if ((current = malloc (sizeof (st_brfile *) * numbrfiles)) == NULL) {
    fprintf (stderr, "malloc current failed, %s\n", strerror (errno));
    exit (-1);
  }

  for (fileno = 0; fileno < numbrfiles; fileno++) {
    if (fileno == 0) {
      get_next_marker (&brfiles[fileno], &data_0);
      do_dkelvin_first_pass (&brfiles[fileno], &data_0);
    } else {
      get_next_marker (&brfiles[fileno], &data_n);
      do_dkelvin_first_pass (&brfiles[fileno], &data_n);
    }
    if (verbose >= 2)
      printf ("first marker in %s is %s\n", brfiles[fileno].name, brfiles[fileno].curmarker.name2);
  }
  print_twopoint_headers (brfiles[0].no_ld);

  while (1) {
    if (mapinfile != NULL) {
      if ((mapptr = next_mapmarker (NULL)) != NULL) {
	strcpy (next_marker.chr, mapptr->chr);
	strcpy (next_marker.name2, mapptr->name);
      } else {
	break;
      }
    } else {
      memcpy (&next_marker, &brfiles[0].curmarker, sizeof (st_brmarker));
    }

    numcurrent = 0;
    alldone = 1;
    for (fileno = 0; fileno < numbrfiles; fileno++) {
      if (brfiles[fileno].eof) {
	if (verbose >= 3)
	  printf ("%s is at EOF\n", brfiles[fileno].name);
	continue;
      }
      alldone = 0;
      if (compare_markers (&next_marker, &brfiles[fileno].curmarker) == -1) {
	if (mapinfile != NULL)
	  continue;
	fprintf (stderr, "marker mismatch at line %d in '%s', expecting %s, found %s\n",
		 brfiles[fileno].lineno, brfiles[fileno].name, next_marker.name2,
		 brfiles[fileno].curmarker.name2);
	exit (-1);
      }
      current[numcurrent++] = &brfiles[fileno];
      compare_headers (current[0], &brfiles[fileno]);
    }
    if (alldone)
      break;
    if (numcurrent == 0) {
      fprintf (stderr, "warning: marker %s from %s doesn't appear in any BR files, skipping\n",
	       next_marker.name2, mapinfile);
      continue;
    }
    if (verbose >= 1) {
      printf ("current marker is '%s', using BR file%s %s", next_marker.name2,
	      (numcurrent > 1) ? "s" : "", current[0]->name);
      for (fileno = 1; fileno < numcurrent; fileno++) 
	printf (", %s", current[fileno]->name);
      printf ("\n");
    }
    
    if (partout != NULL)
      print_partial_header (current[0]);

    sampleno = 0;
    memset (&ldval, 0, sizeof (st_ldvals));
    while (get_data_line (current[0], &data_0) == 1) {
      compare_samples (&data_0, sampleno, current[0]);
      for (fileno = 1; fileno < numcurrent; fileno++) {
	if ((ret = get_data_line (current[fileno], &data_n)) != 1) {
	  fprintf (stderr, "expected data in '%s' at line %d, found %s\n",
		   current[fileno]->name, current[fileno]->lineno,
		   (ret == 0) ? "end-of-file" : "marker line");
	  exit (-1);
	}
	compare_samples (&data_n, sampleno, current[fileno]);
	data_0.lr *= data_n.lr;
      }
      if (partout != NULL)
	print_partial_data (current[0], &data_0);
      
      if (sampleno < 5) {
        ldval.le_small_theta += data_0.lr * dcuhre2[sampleno][2];
      } else if (sampleno < 10) {
        ldval.le_big_theta += data_0.lr * dcuhre2[sampleno][2];
      } else if (sampleno < 140){
        if (dcuhre2[sampleno][1] < cutoff) {
          ldval.ld_small_theta += data_0.lr * dcuhre2[sampleno][2];
        } else {
          ldval.ld_big_theta += data_0.lr * dcuhre2[sampleno][2];
        }
      } else {
        ldval.le_unlinked += data_0.lr * dcuhre2[sampleno][2];
      }
      sampleno++;
    }
    ldval.le_small_theta *= weight;
    ldval.ld_small_theta *= weight;
    ldval.le_big_theta *= (1 - weight);
    ldval.ld_big_theta *= (1 - weight);
    if (forcemap) {
      strcpy (current[0]->curmarker.chr, mapptr->chr);
      current[0]->curmarker.pos = mapptr->avgpos;
    }
    print_twopoint_stats (current[0]->no_ld, &(current[0]->curmarker), &ldval);

    for (fileno = 0; fileno < numcurrent; fileno++) {
      if (&brfiles[0] == current[fileno]) {
	get_next_marker (current[fileno], &data_0);
      } else {
	get_next_marker (current[fileno], &data_n);
      }
    }
  }

  free (current);
  if (data_0.dprimesize > 0)
    free (data_0.dprimes);
  if (data_0.thetasize > 0)
    free (data_0.thetas);
  if (data_n.dprimesize > 0)
    free (data_n.dprimes);
  if (data_n.thetasize > 0)
    free (data_n.thetas);
  for (fileno = 0; fileno < numbrfiles; fileno++) {
    if (brfiles[fileno].datacolsize > 0)
      free (brfiles[fileno].datacols);
  }
  return;
}


void print_twopoint_headers (int no_ld)
{
  if (no_ld) {
    printf ("Chr Trait Marker Position PPL\n");
  } else {
    if (pplinfile == NULL) {
      printf ("Chr Trait Marker Position PPL PPL(LD) PPLD|L PPLD(L)\n");
      if (bfout != NULL)
	fprintf (bfout, "Chr Trait Marker Position PPL(LD) PPLD|L PPLD(L)\n");
    } else if (! allstats) {
      printf ("Chr Trait Marker Position PPLD(L) iPPL iPPLD(L)\n");
      if (bfout != NULL)
	fprintf (bfout, "Chr Trait Marker Position PPLD(L) iPPLD(L)\n");
    } else {
      printf ("Chr Trait Marker Position PPL PPL(LD) PPLD|L PPLD(L) iPPL iPPLD(L)\n");
      if (bfout != NULL)
	fprintf (bfout, "Chr Trait Marker Position PPL(LD) PPLD|L PPLD(L) iPPLD(L)\n");
    }
  }
  if (sixout != NULL)
    fprintf (sixout, "Chr Trait Marker Position LDSmallTheta LDBigTheta LDUnlinked LESmallTheta LEBigTheta LEUnlinked\n");
  return;
}


void print_twopoint_stats (int no_ld, st_brmarker *marker, st_ldvals *ldval)
{
  double ldprior, ldstat;

  printf ("%s %s %s %.4f", marker->chr, marker->name1, marker->name2, marker->pos);

  if (no_ld) {
    printf (" %.3f", calc_upd_ppl (ldval));
  } else {
    if (bfout != NULL)
      fprintf (bfout, "%s %s %s %.4f", marker->chr, marker->name1, marker->name2, marker->pos);
    
    if ((pplinfile == NULL) || (allstats)) {
      printf (" %.3f", calc_upd_ppl (ldval));
      ldstat = calc_upd_ppl_allowing_ld (ldval, DEFAULT_LDPRIOR);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_upd_ppld_given_linkage (ldval, DEFAULT_LDPRIOR);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
    }
    ldstat = calc_upd_ppld_allowing_l (ldval, DEFAULT_LDPRIOR);
    printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
    if (pplinfile != NULL) {
      if ((ldprior = get_ippl (marker->chr, marker->pos)) < MIN_PRIOR)
	ldprior = MIN_PRIOR;
      if (verbose >= 2)
	printf ("ippl is %.6e\n", ldprior);
      ldstat = calc_upd_ppld_allowing_l (ldval, ldprior);
      printf (" %.4f %.*f", ldprior, ldstat >= .025 ? 2 : 4, KROUND (ldstat));
    }
    
    if (bfout != NULL)
      fprintf (bfout, "\n");
  }
  printf ("\n");
  if (sixout != NULL)
    fprintf (sixout, "%s %d %s %.4f %.6e %.6e %.6e %.6e %.6e %.6e\n",  marker->chr,
	     marker->num, marker->name2, marker->pos, ldval->ld_small_theta,
	     ldval->ld_big_theta, ldval->ld_unlinked, ldval->le_small_theta,
	     ldval->le_big_theta, ldval->le_unlinked);
  return;
}


void do_first_pass (st_brfile *brfile, st_multidim *dprimes, st_multidim *thetas, st_data *data)
{
  int ret, lineno;
  long startofdata;

  lineno = brfile->lineno;
  if ((startofdata = ftell (brfile->fp)) == -1) {
    fprintf (stderr, "ftell on file '%s' failed, %s\n", brfile->name, strerror (errno));
    exit (-1);
  }

  while ((ret = get_data_line (brfile, data)) == 1) {
    if (multi_insert (dprimes, data->dprimes, brfile->numdprimes) == -1) {
      fprintf (stderr, "insert into dprimes failed, %s\n", strerror (errno));
      exit (-1);
    }
    if (multi_insert (thetas, data->thetas, brfile->numthetas) == -1) {
      fprintf (stderr, "insert into thetas failed, %s\n", strerror (errno));
      exit (-1);
    }
  }
  if (ret == -1) {
    fprintf (stderr, "can't parse data, line %d in file '%s'\n", brfile->lineno, brfile->name);
    exit (-1);
  }
  
  if (fseek (brfile->fp, startofdata, SEEK_SET) == -1) {
    fprintf (stderr, "fseek on file '%s' failed, %s\n", brfile->name, strerror (errno));
    exit (-1);
  }
  if (brfile->eof)
    brfile->eof = 0;
  brfile->lineno = lineno;
  return;
}


 void do_dkelvin_first_pass (st_brfile *brfile, st_data *data)
 {
   int firstdata, lineno, sampleno, ret;

   if ((firstdata = ftell (brfile->fp)) == -1) {
     fprintf (stderr, "ftell on file '%s' failed, %s\n", brfile->name, strerror (errno));
     exit (-1);
   }
   lineno = brfile->lineno;
   
   sampleno = 0;
   while ((ret = get_data_line (brfile, data)) == 1) {
     compare_samples (data, sampleno, brfile);
     sampleno++;
   }
   
   if (sampleno == 141) 
     brfile->no_ld = 0;
   else if (sampleno != 10) {
     fprintf (stderr, "unexpected number of samples %d in file '%s'\n", sampleno, brfile->name);
     exit (-1);
   } else {
     brfile->no_ld = 1;
   }
   
   if (fseek (brfile->fp, firstdata, SEEK_SET) == -1) {
     fprintf (stderr, "fseek on file '%s' failed, %s\n", brfile->name, strerror (errno));
     exit (-1);
   }
   if (brfile->eof)
     brfile->eof = 0;
   brfile->lineno = lineno;
   return;
 }



void calc_ldvals_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr, st_ldvals *ldval)
{
  int didx, zero_didx, thidx, ndprimes = 0, va;
  double mtheta1, mtheta2, lr1, lr2, cutlr, *zeros;
  st_dim *mthetas;
 
  ldval->ld_small_theta = ldval->ld_big_theta = ldval->ld_unlinked = 0;
  ldval->le_small_theta = ldval->le_big_theta = ldval->le_unlinked = 0;

  if ((zeros = malloc (sizeof (double) * dprimes->numdims)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes->numdims; va++)
    zeros[va] = 0.0;
  if ((zero_didx = multi_find (dprimes, zeros, dprimes->numdims)) == -1) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }
  free (zeros);
  mthetas = &thetas->dims[0];

  for (didx = 0; didx < dprimes->totalelems; didx++) {

    /* See comment in calc_ppl_sexavg() */
    if (lr[didx][0] == DBL_MAX)
      continue;
    
    if (didx == zero_didx) {
      /* D' == 0 */
      for (thidx = 1; thidx < mthetas->numelems; thidx++) {
	mtheta1 = mthetas->arr[thidx - 1];
	mtheta2 = mthetas->arr[thidx];
        lr1 = lr[didx][thidx-1];
        lr2 = lr[didx][thidx];
	
	if (mtheta1 >= cutoff) {
	  ldval->le_big_theta += (lr1 + lr2) * (mtheta2 - mtheta1);
	} else if (mtheta2 <= cutoff) {
	  ldval->le_small_theta += (lr1 + lr2) * (mtheta2 - mtheta1);
	} else {
	  cutlr = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	  ldval->le_small_theta += (lr1 + cutlr) * (cutoff - mtheta1);
	  ldval->le_big_theta += (cutlr + lr2) * (mtheta2 - cutoff);
	}
	if (mtheta2 == 0.5)
	  ldval->le_unlinked += lr2;
      }
    } else {
      /* D' != 0 */
      ndprimes++;
      for (thidx = 1; thidx < mthetas->numelems; thidx++) {	
	mtheta1 = mthetas->arr[thidx - 1];
	mtheta2 = mthetas->arr[thidx];
        lr1 = lr[didx][thidx-1];
        lr2 = lr[didx][thidx];
	
	if (mtheta1 >= cutoff) {
	  ldval->ld_big_theta += (lr1 + lr2) * (mtheta2 - mtheta1);
	} else if (mtheta2 <= cutoff) {
	  ldval->ld_small_theta += (lr1 + lr2) * (mtheta2 - mtheta1);
	} else {
	  cutlr = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	  ldval->ld_small_theta += (lr1 + cutlr) * (cutoff - mtheta1);
	  ldval->ld_big_theta += (cutlr + lr2) * (mtheta2 - cutoff);
	}
	if (mtheta2 == 0.5)
	  ldval->ld_unlinked += lr2;
      }
    }
  }

  /* The 0.5 here is factored out of the calculation of area */
  ldval->ld_small_theta *= 0.5 * (weight / cutoff);
  ldval->le_small_theta *= 0.5 * (weight / cutoff);
  ldval->ld_big_theta *= 0.5 * ((1 - weight) / (0.5 - cutoff));
  ldval->le_big_theta *= 0.5 * ((1 - weight) / (0.5 - cutoff));

  /* Average LD values over the number of non-zero D's */
  ldval->ld_small_theta /= ndprimes;
  ldval->ld_big_theta /= ndprimes;
  ldval->ld_unlinked /= ndprimes;

  return;
}


void calc_ldvals_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr, st_ldvals *ldval)
{
  int didx, mthidx, fthidx, va, ndprimes = 0, zero_didx;
  double mtheta1, mtheta2, ftheta1, ftheta2, lr1, lr2, lr3, lr4, cutlr2, cutlr3, cutlr4;
  double cutvol, *zeros;
  st_dim *mthetas, *fthetas;
 
  ldval->ld_small_theta = ldval->ld_big_theta = ldval->ld_unlinked = 0;
  ldval->le_small_theta = ldval->le_big_theta = ldval->le_unlinked = 0;

  if ((zeros = malloc (sizeof (double) * dprimes->numdims)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes->numdims; va++)
    zeros[va] = 0.0;
  if ((zero_didx = multi_find (dprimes, zeros, dprimes->numdims)) == -1) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }
  free (zeros);
  mthetas = &thetas->dims[0];
  fthetas = &thetas->dims[1];

  for (didx = 0; didx < dprimes->totalelems; didx++) {

    /* See comment in calc_ppl_sexavg() */
    if (lr[didx][0] == DBL_MAX)
      continue;

    if (didx == zero_didx) {
      /* D' == 0 */
      for (mthidx = 1; mthidx < mthetas->numelems; mthidx++) {
	mtheta1 = mthetas->arr[mthidx - 1];
	mtheta2 = mthetas->arr[mthidx];
	for (fthidx = 1; fthidx < fthetas->numelems; fthidx++) {
	  ftheta1 = fthetas->arr[fthidx - 1];
	  ftheta2 = fthetas->arr[fthidx];
	  lr1 = lr[didx][(mthidx - 1) * mthetas->numelems + fthidx - 1];
	  lr2 = lr[didx][mthidx * mthetas->numelems + fthidx - 1];
	  lr3 = lr[didx][(mthidx - 1) * mthetas->numelems + fthidx];
	  lr4 = lr[didx][mthidx * mthetas->numelems + fthidx];
	  
	  if ((mtheta1 >= cutoff) || (ftheta1 >= cutoff)) {
	    /* entire region is outside cutoff area */
	    ldval->le_big_theta += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) *
	      (lr1 + lr2 + lr3 + lr4);
	    
	  } else if ((mtheta2 <= cutoff) && (ftheta2 <= cutoff)) {
	    /* entire region inside cutoff area */
	    ldval->le_small_theta += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) *
	      (lr1 + lr2 + lr3 + lr4);
	    
	  } else if ((ftheta2 <= cutoff) && ((mtheta1 <= cutoff) && (cutoff <= mtheta2))) {
	    /* region straddles cutoff on male axis */
	    cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	    cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	    ldval->le_small_theta += (cutoff - mtheta1) * (ftheta2 - ftheta1) *
	      (lr1 + cutlr2 + lr3 + cutlr4);
	    ldval->le_big_theta += (mtheta2 - cutoff) * (ftheta2 - ftheta1) *
	      (cutlr2 + lr2 + cutlr4 + lr4);
	    
	  } else if ((mtheta2 <= cutoff) && ((ftheta1 <= cutoff) && (cutoff <= ftheta2))) {
	    /* region straddles cutoff on female axis */
	    cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	    cutlr4 = lr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr4 - lr2);
	    ldval->le_small_theta += (cutoff - ftheta1) * (mtheta2 - mtheta1) *
	      (lr1 + lr2 + cutlr3 + cutlr4);
	    ldval->le_big_theta += (ftheta2 - cutoff) * (mtheta2 - mtheta1) *
	      (cutlr3 + cutlr4 + lr3 + lr4);
	    
	  } else {
	    /* region straddles cutoff on both axes */
	    ldval->le_big_theta += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) *
	      (lr1 + lr2 + lr3 + lr4);
	    cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	    cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	    cutlr4 = cutlr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (cutlr4 - cutlr2);
	    cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	    cutvol = (cutoff - mtheta1) * (cutoff - ftheta1) * (lr1 + cutlr2 + cutlr3 + cutlr4);
	    ldval->le_big_theta -= cutvol;
	    ldval->le_small_theta += cutvol;
	  }
	  if ((mtheta2 == 0.5) && (ftheta2 == 0.5))
	    ldval->le_unlinked += lr4;
	}
      }
    } else {
      /* D' != 0 */
      ndprimes++;
      for (mthidx = 1; mthidx < mthetas->numelems; mthidx++) {
	mtheta1 = mthetas->arr[mthidx - 1];
	mtheta2 = mthetas->arr[mthidx];
	for (fthidx = 1; fthidx < fthetas->numelems; fthidx++) {
	  ftheta1 = fthetas->arr[fthidx - 1];
	  ftheta2 = fthetas->arr[fthidx];
	  lr1 = lr[didx][(mthidx - 1) * mthetas->numelems + fthidx - 1];
	  lr2 = lr[didx][mthidx * mthetas->numelems + fthidx - 1];
	  lr3 = lr[didx][(mthidx - 1) * mthetas->numelems + fthidx];
	  lr4 = lr[didx][mthidx * mthetas->numelems + fthidx];
	  
	  if ((mtheta1 >= cutoff) || (ftheta1 >= cutoff)) {
	    /* entire region is outside cutoff area */
	    ldval->ld_big_theta += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) *
	      (lr1 + lr2 + lr3 + lr4);
	    
	  } else if ((mtheta2 <= cutoff) && (ftheta2 <= cutoff)) {
	    /* entire region inside cutoff area */
	    ldval->ld_small_theta += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) *
	      (lr1 + lr2 + lr3 + lr4);
	    
	  } else if ((ftheta2 <= cutoff) && ((mtheta1 <= cutoff) && (cutoff <= mtheta2))) {
	    /* region straddles cutoff on male axis */
	    cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	    cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	    ldval->ld_small_theta += (cutoff - mtheta1) * (ftheta2 - ftheta1) *
	      (lr1 + cutlr2 + lr3 + cutlr4);
	    ldval->ld_big_theta += (mtheta2 - cutoff) * (ftheta2 - ftheta1) *
	      (cutlr2 + lr2 + cutlr4 + lr4);
	    
	  } else if ((mtheta2 <= cutoff) && ((ftheta1 <= cutoff) && (cutoff <= ftheta2))) {
	    /* region straddles cutoff on female axis */
	    cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	    cutlr4 = lr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr4 - lr2);
	    ldval->ld_small_theta += (cutoff - ftheta1) * (mtheta2 - mtheta1) *
	      (lr1 + lr2 + cutlr3 + cutlr4);
	    ldval->ld_big_theta += (ftheta2 - cutoff) * (mtheta2 - mtheta1) *
	      (cutlr3 + cutlr4 + lr3 + lr4);
	    
	  } else {
	    /* region straddles cutoff on both axes */
	    ldval->ld_big_theta += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) *
	      (lr1 + lr2 + lr3 + lr4);
	    cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	    cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	    cutlr4 = cutlr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (cutlr4 - cutlr2);
	    cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	    cutvol = (cutoff - mtheta1) * (cutoff - ftheta1) * (lr1 + cutlr2 + cutlr3 + cutlr4);
	    ldval->ld_big_theta -= cutvol;
	    ldval->ld_small_theta += cutvol;
	  }
	  if ((mtheta2 == 0.5) && (ftheta2 == 0.5))
	    ldval->ld_unlinked += lr4;
	}
      }
    }
  }

  /* The first 0.25 here is factored out of the calculation of volume */
  /* The second 0.25 is the square of the largest value of Theta (0.5) */
  ldval->ld_small_theta *= 0.25 * (weight / (cutoff * cutoff));
  ldval->le_small_theta *= 0.25 * (weight / (cutoff * cutoff));
  ldval->ld_big_theta *= 0.25 * ((1 - weight) / (0.25 - (cutoff * cutoff)));
  ldval->le_big_theta *= 0.25 * ((1 - weight) / (0.25 - (cutoff * cutoff)));

  /* Average LD values over the number of non-zero D's */
  ldval->ld_small_theta /= ndprimes;
  ldval->ld_big_theta /= ndprimes;
  ldval->ld_unlinked /= ndprimes;

  return;
}


 double calc_upd_ppl (st_ldvals *ldval)
 {
   double numerator;
   double denomRight;
   double ppl;

   numerator = 
     ldval->le_small_theta * prior +
     ldval->le_big_theta * prior;
   denomRight = 1 - prior;
   ppl = numerator / (numerator + denomRight);
 }


double calc_upd_ppl_allowing_ld (st_ldvals *ldval, double ldprior)
{
  double numerator;
  double denomRight;
  double ldppl;

  numerator =
    ldval->ld_small_theta * ldprior * 0.021 + 
    ldval->ld_big_theta * ldprior * 0.0011+ 
    ldval->le_small_theta * ldprior * 0.979 + 
    ldval->le_big_theta * ldprior * 0.9989;
  denomRight =
    ldval->le_unlinked * (1 - ldprior);
  ldppl = numerator / (numerator + denomRight);
  
  if (bfout != NULL) 
    fprintf (bfout, " %.6e", log10 (numerator / denomRight));
  return (ldppl);
}


double calc_upd_ppld_given_linkage (st_ldvals *ldval, double ldprior)
{
  double numerator;
  double denomRight;
  double ppld_given_l;

  numerator =
    ldval->ld_small_theta * ldprior * 0.021 + 
    ldval->ld_big_theta * ldprior * 0.0011;
  denomRight =
    ldval->le_small_theta * ldprior * 0.979 + 
    ldval->le_big_theta * ldprior * 0.9989;
  ppld_given_l = numerator / (numerator + denomRight);

  if (bfout != NULL) 
    fprintf (bfout, " %.6e", log10 (numerator / denomRight));
  return (ppld_given_l);
}


double calc_upd_ppld_allowing_l (st_ldvals *ldval, double ldprior)
{
  double numerator;
  double denomRight;
  double ppld; 

  numerator =
    ldval->ld_small_theta * ldprior * 0.021 +
    ldval->ld_big_theta * ldprior * 0.0011; 
  denomRight =
    ldval->le_small_theta * ldprior * 0.979 + 
    ldval->le_big_theta * ldprior * 0.9989 + 
    ldval->le_unlinked * (1 - ldprior);
  ppld = numerator / (numerator + denomRight); 

  if (bfout != NULL) 
    fprintf (bfout, " %.6e", log10 (numerator / denomRight));
  return (ppld);
}


#define OPT_SEXSPEC   1
#define OPT_MULTI     2
#define OPT_RELAX     3
#define OPT_ALLSTATS  4
#define OPT_OKELVIN   5
#define OPT_VERBOSE   6
#define OPT_PRIOR     7
#define OPT_WEIGHT    8
#define OPT_CUTOFF    9
#define OPT_MAPIN    10
#define OPT_PPLIN    11
#define OPT_PARTOUT  12
#define OPT_BFOUT    13
#define OPT_SIXOUT   14
#define OPT_FORCEMAP 15
#define OPT_HELP     16

int parse_command_line (int argc, char **argv)
{
  int arg, long_arg, long_idx;
  char *partoutfile=NULL, *bfoutfile=NULL, *sixoutfile=NULL;
  struct option cmdline[] = { { "sexspecific", 0, &long_arg, OPT_SEXSPEC },
			      { "multipoint", 0, &long_arg, OPT_MULTI },
			      { "relax", 1, &long_arg, OPT_RELAX },
			      { "allstats", 0, &long_arg, OPT_ALLSTATS },
			      { "okelvin", 0, &long_arg, OPT_OKELVIN },
			      { "verbose", 1, &long_arg, OPT_VERBOSE },
			      { "prior", 1, &long_arg, OPT_PRIOR },
			      { "weight", 1, &long_arg, OPT_WEIGHT },
			      { "cutoff", 1, &long_arg, OPT_CUTOFF },
			      { "pplin", 1, &long_arg, OPT_PPLIN },
			      { "mapin", 1, &long_arg, OPT_MAPIN },
			      { "partout", 1, &long_arg, OPT_PARTOUT },
			      { "bfout", 1, &long_arg, OPT_BFOUT },
			      { "sixout", 1, &long_arg, OPT_SIXOUT },
			      { "forcemap", 0, &long_arg, OPT_FORCEMAP },
			      { "help", 0, &long_arg, OPT_HELP },
			      { NULL, 0, NULL, 0 } };
  struct stat statbuf;
  
  while ((arg = getopt_long (argc, argv, "smraovp:w:c:", cmdline, &long_idx)) != -1) {
    if ((arg == 's') || ((arg == 0) && (long_arg == OPT_SEXSPEC))) {
      sexspecific = 1;

    } else if ((arg == 'm') || ((arg == 0) && (long_arg == OPT_MULTI))) {
      multipoint = 1;

    } else if ((arg == 'r') || ((arg == 0) && (long_arg == OPT_RELAX))) {
      relax = 1;

    } else if ((arg == 'a') || ((arg == 0) && (long_arg == OPT_ALLSTATS))) {
      allstats = 1;

    } else if ((arg == 'o') || ((arg == 0) && (long_arg == OPT_OKELVIN))) {
      okelvin = 1;

    } else if ((arg == 'v') || ((arg == 0) && (long_arg == OPT_VERBOSE))) {
      verbose++;

    } else if (arg == 'p') {
      prior = validate_double_arg (optarg, "-p");
    } else if ((arg == 0) && (long_arg == OPT_PRIOR)) {
      prior = validate_double_arg (optarg, "--prior");

    } else if (arg == 'w') {
      weight = validate_double_arg (optarg, "-w");
    } else if ((arg == 0) && (long_arg == OPT_WEIGHT)) {
      weight = validate_double_arg (optarg, "--weight");

    } else if (arg == 'c') {
      cutoff = validate_double_arg (optarg, "-c");
    } else if ((arg == 0) && (long_arg == OPT_CUTOFF)) {
      cutoff = validate_double_arg (optarg, "--cutoff");

    } else if ((arg == 0) && (long_arg == OPT_MAPIN)) {
      mapinfile = optarg;

    } else if ((arg == 0) && (long_arg == OPT_PPLIN)) {
      pplinfile = optarg;

    } else if ((arg == 0) && (long_arg == OPT_PARTOUT)) {
      partoutfile = optarg;

    } else if ((arg == 0) && (long_arg == OPT_BFOUT)) {
      bfoutfile = optarg;

    } else if ((arg == 0) && (long_arg == OPT_SIXOUT)) {
      sixoutfile = optarg;

    } else if ((arg == 0) && (long_arg == OPT_FORCEMAP)) {
      forcemap = 1;

    } else if ((arg == 0) && (long_arg == OPT_HELP)) {
      usage ();
      exit (0);

    } else {
      /* getopt_long will emit a message describing the bad option/argument */
      usage ();
      exit (-1);
    }
  }

  if (prior < MIN_PRIOR) {
    fprintf (stderr, "warning: specified prior too small, fixing prior at %.4e\n", MIN_PRIOR);
    prior = MIN_PRIOR;
  }

  if ((multipoint) && (bfoutfile != NULL)) {
    fprintf (stderr, "%s: --bfout is nonsensical with --multipoint\n", pname);
    exit (-1);
  }

  if ((multipoint) && (partoutfile != NULL)) {
    fprintf (stderr, "%s: --partout is nonsensical with --multipoint\n", pname);
    exit (-1);
  }

  if ((multipoint) && (pplinfile != NULL)) {
    fprintf (stderr, "%s: --pplin is nonsensical with --multipoint\n", pname);
    exit (-1);
  }

  if ((multipoint) && (allstats)) {
    fprintf (stderr, "%s: --allstats is nonsensical with --multipoint\n", pname);
    exit (-1);
  }

  if ((mapinfile != NULL) && (relax)) {
    fprintf (stderr, "%s: --relax cannot be combined with --mapin\n", pname);
    exit (-1);
  }

  if ((forcemap) && (mapinfile == NULL)) {
    fprintf (stderr, "%s: --forcemap requires --mapin\n", pname);
    exit (-1);
  }

  if (partoutfile != NULL) {
    if (stat (partoutfile, &statbuf) != -1) {
      fprintf (stderr, "%s: won't open '%s' for writing, file exists\n", pname, partoutfile);
      exit (-1);
    }
    if ((partout = fopen (partoutfile, "w")) == NULL) {
      fprintf (stderr, "%s: open '%s' for writing failed, %s\n", pname, partoutfile,
	       strerror (errno));
      exit (-1);
    }
  }

  if (bfoutfile != NULL) {
    /*if (stat (bfoutfile, &statbuf) != -1) {
      fprintf (stderr, "%s: won't open '%s' for writing, file exists\n", pname, bfoutfile);
      exit (-1);
      }*/
    if ((bfout = fopen (bfoutfile, "a")) == NULL) {
      fprintf (stderr, "%s: open '%s' for writing failed, %s\n", pname, bfoutfile,
	       strerror (errno));
      exit (-1);
    }
  }

  if (sixoutfile != NULL) {
    if (stat (sixoutfile, &statbuf) != -1) {
      fprintf (stderr, "%s: won't open '%s' for writing, file exists\n", pname, sixoutfile);
      exit (-1);
    }
    if ((sixout = fopen (sixoutfile, "w")) == NULL) {
      fprintf (stderr, "%s: open '%s' for writing failed, %s\n", pname, sixoutfile,
	       strerror (errno));
      exit (-1);
    }
  }
  
  return (optind);
}


void usage ()
{
  printf ("usage: %s [ options ] brfile [brfile...]\n  valid options:\n", pname);
  printf ("  -s|--sexspecific : input data contains sex-specific Thetas\n");
  printf ("  -m|--multipoint : input data is multipoint\n");
  printf ("  -o|--okelvin : input data is from original (fixed-grid) kelvin\n");
  printf ("  -r|--relax : supress comparing marker names across input files\n");
  printf ("  -p <num>|--prior <num> : set linkage prior probability to <num>\n");
  printf ("  -c <num>|--cutoff <num> : set small-Theta cutoff to <num>\n");
  printf ("  -w <num>|--wieght <num> : set small-Theta weight to <num>\n");
  printf ("  --mapin <mapfile> : use mapfile to order markers\n");
  printf ("  --partout <partfile> : write updated Bayes Ratios to partfile\n");
  printf ("  --help : display this help text\n");
  return;
}


double validate_double_arg (char *arg, char *optname)
{
  double da;
  char *ptr;

  if (((da = strtod (arg, &ptr)) == 0) && (arg == ptr)) {
    fprintf (stderr, "%s: argument '%s' to option '%s' is illegal\n", pname, arg, optname);
    exit (-1);
  }
  return (da);
}


void open_brfile (st_brfile *brfile)
{
  char buff[BUFFLEN];
  int major, minor, patch, fileverno, curverno;
  
  if ((brfile->fp = fopen (brfile->name, "r")) == NULL) {
    fprintf (stderr, "open '%s' failed, %s\n", brfile->name, strerror (errno));
    exit (-1);
  }
  if (fgets (buff, BUFFLEN, brfile->fp) == NULL) {
    if (brfile->eof = feof (brfile->fp))
      fprintf (stderr, "error reading %s at line 0, %s\n", brfile->name, strerror (errno));
    else 
      fprintf (stderr, "file %s is empty\n", brfile->name);
    exit (-1);
  }
  sscanf (current, "%d.%d.%d", &major, &minor, &patch);
  curverno = major * 100000 + minor * 1000 + patch;

  if (sscanf (buff, "# Version V%d.%d.%d", &major, &minor, &patch) != 3) {
    fprintf (stderr, "file %s is unversioned, please convert to at least V%s\n",
	     brfile->name, current);
    exit (-1);
  }
  fileverno = major * 100000 + minor * 1000 + patch;
  if (fileverno < curverno) {
    fprintf (stderr, "file %s is version V%d.%d.%d, please convert to at least V%s\n",
	     brfile->name, major, minor, patch, current);
    exit (-1);
  }
  brfile->lineno = 1;
  return;
}


void get_next_marker (st_brfile *brfile, st_data *data)
{
  int ret, lineno;
  long start_of_data;

  if (get_marker_line (brfile) == 0) {
    memset (&brfile->curmarker, 0, sizeof (st_brmarker));
    return;
  }
  get_header_line (brfile);
  
  if (data->dprimesize < brfile->numdprimes) {
    data->dprimesize = brfile->numdprimes;
    if ((data->dprimes = realloc (data->dprimes, sizeof (double) * data->dprimesize)) == NULL) {
      fprintf (stderr, "realloc failed, %s\n", strerror (errno));
      exit (-1);
    }
  }
  if (data->thetasize < brfile->numthetas) {
    data->thetasize = brfile->numthetas;
    if ((data->thetas = realloc (data->thetas, sizeof (double) * data->thetasize)) == NULL) {
      fprintf (stderr, "realloc failed, %s\n", strerror (errno));
      exit (-1);
    }
  }
  if ((start_of_data = ftell (brfile->fp)) == -1) {
    fprintf (stderr, "ftell failed for '%s', %s\n", brfile->name, strerror (errno));
    exit (-1);
  }
  lineno = brfile->lineno;
  if ((ret = get_data_line (brfile, data)) == 0) {
    fprintf (stderr, "file '%s' ends unexpectedly at line %d\n", brfile->name, brfile->lineno);
    exit (-1);
  } else if (ret == 2) {
    fprintf (stderr, "expected data, found marker at line %d in '%s'\n", brfile->lineno
	     , brfile->name);
    exit (-1);
  }
  if (fseek (brfile->fp, start_of_data, SEEK_SET) == -1) {
    fprintf (stderr, "fseek failed for '%s', %s\n", brfile->name, strerror (errno));
    exit (-1);
  }
  brfile->lineno = lineno;
  return;
}


int get_marker_line (st_brfile *brfile)
{
  char buff[BUFFLEN];
  st_brmarker *marker;

  marker = &brfile->curmarker;
  
  if (fgets (buff, BUFFLEN, brfile->fp) == NULL) {
    /* An end-of-file looking for a marker is not necessarily fatal, so let the caller decide */
    if (brfile->eof = feof (brfile->fp))
      return (0);
    fprintf (stderr, "error reading '%s' at line %d, %s\n", brfile->name, brfile->lineno,
	     strerror (errno));
    exit (-1);
  }
  brfile->lineno++;
  
  if (sscanf (buff, "# %d %s %s", &marker->num, marker->name1, marker->name2) != 3) {
    if (strncasecmp (buff, "Chr", 3) == 0) {
      fprintf (stderr, "expected marker line, line %d in '%s', found possible header line\n",
	       brfile->lineno, brfile->name);
      fprintf (stderr, "maybe %s contains multipoint data?\n", brfile->name);
    } else {
      fprintf (stderr, "can't parse marker, line %d in file '%s'\n", brfile->lineno, brfile->name);
    }
    exit (-1);
  }
  
  return (1);
}


int get_header_line (st_brfile *brfile)
{ 
  char buff[BUFFLEN], token[32];
  char *pa, *pb, *pc=NULL;
  int actualcols, numlrcols=0, va;

  if (fgets (buff, BUFFLEN, brfile->fp) == NULL) {
    if (feof (brfile->fp))
      fprintf (stderr, "file '%s' ends where header line expected at line %d\n",
	       brfile->name, brfile->lineno);
    else 
      fprintf (stderr, "error reading '%s' at line %d, %s\n", brfile->name, brfile->lineno,
	       strerror (errno));
    exit (-1);
  }
  brfile->lineno++;
  brfile->numcols = brfile->numdprimes = brfile->numthetas = brfile->no_ld = 
    brfile->holey_grid = brfile->two_point = 0;

  pa = strtok_r (buff, " \t\n", &pb);
  while (pa != NULL) {
    strcpy (token, pa);
    while ((strchr (token, '(') != NULL) && (strrchr (token, ')') == NULL)) {
      if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL) {
	fprintf (stderr, "can't parse header, line %d in file '%s'\n", brfile->lineno,
		 brfile->name);
	exit (-1);
      }
      strcat (token, pa);
    }
    
    if (token[0] == '(') {
      fprintf (stderr, "header '%s' has no identifying string before '('\n", token);
      exit (-1);
    }
    
    if ((pa = strchr (token, '(')) == NULL) {
      actualcols = 1;
    } else {
      /* Some parenthesized expression */
      actualcols = 0;
      *pa = '\0';
      pa = strtok_r (++pa, ",", &pc);
      while (pa != NULL) {
	actualcols++;
	pa = strtok_r (NULL, ",", &pc);
      }
    }

    brfile->numcols += actualcols;
    if (brfile->numcols > brfile->datacolsize) {
      brfile->datacolsize = brfile->numcols;
      if ((brfile->datacols = realloc (brfile->datacols, sizeof (int) * brfile->datacolsize)) == NULL) {
	fprintf (stderr, "realloc failed, %s\n", strerror (errno));
	exit (-1);
      }
    }
    
    if (((token[0] == 'D') && ((token[1] >= '0') && (token[1] <= '9')) &&
	 ((token[1] >= '0') && (token[1] <= '9'))) ||
	(strcasecmp (token, "DPrime") == 0)) {
      /* A D-prime column */
      brfile->two_point = 1;
      brfile->numdprimes += actualcols;
#ifdef DEBUG
      printf ("number of dprime cols %d\n", brfile->numdprimes);
#endif
      for (; actualcols > 0; actualcols--) 
	brfile->datacols[brfile->numcols - actualcols] = DPRIME_COL;
      
    } else if (strcasecmp (token, "Theta") == 0) {
      /* A Theta column */
      brfile->two_point = 1;
      if (sexspecific) {
	brfile->numthetas += actualcols;
#ifdef DEBUG
	printf ("number of theta cols %d\n", brfile->numthetas);
#endif
	for (; actualcols > 0; actualcols--) 
	  brfile->datacols[brfile->numcols - actualcols] = THETA_COL;
      } else {
	brfile->numthetas++;
#ifdef DEBUG
	printf ("number of theta cols %d\n", brfile->numthetas);
#endif
	brfile->datacols[brfile->numcols - actualcols] = THETA_COL;
	for (--actualcols; actualcols > 0; actualcols--) 
	  brfile->datacols[brfile->numcols - actualcols] = 0;
      }
      
    } else if ((strcasecmp (token, "AVG_LR") == 0) || (strcasecmp (token, "AVGLR") == 0) ||
	       (strcasecmp (token, "BR") == 0) || (strcasecmp (token, "BayesRatio") == 0)) {
      /* The Bayes Ratio column (formerly the average likelihood ratio) */
      brfile->datacols[brfile->numcols - 1] = LR_COL;
      numlrcols++;
#ifdef DEBUG
      printf ("number of lr cols %d\n", numlrcols);
#endif
      
    } else if ((strcasecmp (token, "Pos") == 0) || (strcasecmp (token, "Position") == 0)) {
      /* The position column */
      brfile->datacols[brfile->numcols - 1] = POS_COL;
#ifdef DEBUG
      printf ("position col\n");
#endif
      
    } else if ((strcasecmp (token, "Chr") == 0) || (strcasecmp (token, "Chromosome") == 0)) {
      /* The chromosome column */
      brfile->datacols[brfile->numcols - 1] = CHR_COL;
#ifdef DEBUG
      printf ("chromosome col\n");
#endif
      
    } else if ((strcasecmp (token, "MarkerList") == 0) || (strcasecmp (token, "PPL") == 0)) {
      brfile->two_point = 0;

    } else {
      /* Something else */
      for (; actualcols > 0; actualcols--) 
        brfile->datacols[brfile->numcols - actualcols] = 0;
    }

#ifdef DEBUG
    printf ("\n");
#endif
    pa = strtok_r (NULL, " \t\n", &pb);
  }

  if (numlrcols != 1) {
    fprintf (stderr, "header at line %d in '%s' has no BR column\n", brfile->lineno,
	     brfile->name);
    exit (-1);
  }

  if (brfile->numdprimes == 0) {
    /* If there's no D' columns, then we've got a non-LD run, so we set that flag */
    brfile->no_ld = 1;
    if (brfile->two_point) {
      /* We need at least one dimension of D's for our 2-dimensional array of avgLRs,
       * though, so we'll fake D' == 0 for every data line, and we need allocate
       * space for our dummy D'.
       */
      brfile->numdprimes = 1;
    }
  } else if (brfile->numdprimes > 1) {
    brfile->holey_grid = 1;
  }
  
  #ifdef DEBUG
  for (numlrcols = 0; numlrcols < brfile->numcols; numlrcols++) {
    printf ("col %d is type %d\n", numlrcols, brfile->datacols[numlrcols]);
  }
  #endif

  return (1);
}


int get_data_line (st_brfile *brfile, st_data *data)
{
  char buff[BUFFLEN], *endptr=NULL;
  char *pa=NULL, *pb=NULL;
  long previous;
  int va = 0, dprimecnt=0, thetacnt=0;
  st_brmarker *marker;

  marker = &brfile->curmarker;
  if ((previous = ftell (brfile->fp)) == -1) {
    fprintf (stderr, "ftell failed, %s\n", strerror (errno));
    exit (-1);
  }
  if (fgets (buff, BUFFLEN, brfile->fp) == NULL) {
    if (brfile->eof = feof (brfile->fp))
      return (0);
    fprintf (stderr, "error reading '%s' at line %d, %s\n", brfile->name, brfile->lineno,
	     strerror (errno));
    exit (-1);
  }
  
  if ((pa = strtok_r (buff, " (),\t\n", &pb)) == NULL) {
    fprintf (stderr, "unexpectedly short line %d in '%s'\n", brfile->lineno+1, brfile->name);
    exit (-1);
  }

  /* If it's a new marker line, back the file pointer up, and return a 
   * non-error, non-dataline value.
   */
  if (strcmp (pa, "#") == 0) {
    if (fseek (brfile->fp, previous, SEEK_SET) == -1) {
      fprintf (stderr, "fseek failed, %s\n", strerror (errno));
      exit (-1);
    }
    return (2);
  }
  brfile->lineno++;

  while (1) {
    switch (brfile->datacols[va]) {
    case DPRIME_COL:
      data->dprimes[dprimecnt++] = strtod (pa, &endptr);
      break;
    case THETA_COL:
      data->thetas[thetacnt++] = strtod (pa, &endptr);
      break;
    case LR_COL:
      data->lr = strtod (pa, &endptr);
      break;
    case POS_COL:
      marker->pos = strtod (pa, &endptr);
      break;
    case CHR_COL:
      endptr = strcpy (marker->chr, pa);
      break;
    }
    if (pa == endptr) {
      fprintf (stderr, "illegal data in line %d in '%s'\n", brfile->lineno, brfile->name);
      exit (-1);
    }
    if (++va >= brfile->numcols)
      break;
    if ((pa = strtok_r (NULL, " (),\t\n", &pb)) == NULL) {
      fprintf (stderr, "unexpectedly short line %d in '%s'\n", brfile->lineno, brfile->name);
      exit (-1);
    }
  }

  /* If we've got a two-point, non-LD input file, we need to dummy up a D' of 0
   * for each data line, just so we've got a D' dimension for our
   * 2-dimensional array of avgLRs.
   */
  if ((brfile->two_point) && (brfile->no_ld))
    data->dprimes[0] = 0;

#ifdef DEBUG
  printf ("%5d: chr %s, pos %6.4f, dprimes", brfile->lineno, marker->chr, marker->pos);
  for (va = 0; va < dprimecnt; va++) {
    printf (" %5.2f", data->dprimes[va]);
  }
  printf (", thetas");
  for (va = 0; va < thetacnt; va++) {
    printf (" %5.2f", data->thetas[va]);
  }
  printf (", BR %8.6e, buff %lx, pa %lx, pb %lx\n", data->lr, buff, pa, pb);
#endif

  return (1);
}


void print_partial_header (st_brfile *brfile)
{
  int colno=0, dprimecnt=0;
  st_brmarker *marker;

  marker = &brfile->curmarker;
  fprintf (partout, "# %d %s %s\n", marker->num, marker->name1, marker->name2);
  while (colno < brfile->numcols) {
    switch (brfile->datacols[colno]) {
    case CHR_COL:
      /* It's assumed here that the chromosome column will be first (no leading space) */
      fprintf (partout, "Chr");
      colno++;
      break;
    case POS_COL:
      fprintf (partout, " Position");
      colno++;
      break;
    case DPRIME_COL:
      fprintf (partout, " D1%d", ++dprimecnt);
      colno++;
      break;
    case THETA_COL:
      fprintf (partout, " Theta(M,F)");
      colno += (sexspecific) ? 2 : 1;
      break;
    case LR_COL:
      fprintf (partout, " BayesRatio");
      colno++;
      break;
    default:
      colno++;
    }
  }
  fprintf (partout, "\n");
}


void print_partial_data (st_brfile *brfile, st_data *data)
{
  int colno=0, dprimecnt=0;
  st_brmarker *marker;

  marker = &brfile->curmarker;
  while (colno < brfile->numcols) {
    switch (brfile->datacols[colno]) {
    case CHR_COL:
      fprintf (partout, "%s", marker->chr);
      colno++;
      break;
    case POS_COL:
      fprintf (partout, " %.4f", marker->pos);
      colno++;
      break;
    case DPRIME_COL:
      fprintf (partout, " %.2f", data->dprimes[dprimecnt++]);
      colno++;
      break;
    case THETA_COL:
      if (! sexspecific) {
	fprintf (partout, " (%.4f,%.4f)", data->thetas[0], data->thetas[0]);
	colno++;
      } else {
	fprintf (partout, " (%.4f,%.4f)", data->thetas[0], data->thetas[1]);
	colno += 2;
      }
      break;
    case LR_COL:
      fprintf (partout, " %.6e", data->lr);
      colno++;
      break;
    default:
      colno++;
    }
    
  }
  fprintf (partout, "\n");
}


void compare_headers (st_brfile *f1, st_brfile *f2)
{
  if ((f1->numdprimes != f2->numdprimes) || (f1->no_ld != f2->no_ld)) {
    fprintf (stderr, "header at line %d in '%s' has %d D' columns, expected %d\n",
	     f2->lineno, f2->name, (f2->no_ld == 0) ? f2->numdprimes : 0,
	     (f1->no_ld == 0) ? f1->numdprimes : 0);
    exit (-1);
  }
  if (f1->numthetas != f2->numthetas) {
    fprintf (stderr, "header at line %d in '%s' has %d Theta columns, expected %d\n",
	     f2->lineno, f2->name, f2->numthetas, f1->numthetas);
    exit (-1);
  }
  if (f1->two_point != f2->two_point) {
    fprintf (stderr, "file '%s' contains %s data, expected %s data\n", f2->name,
	     (f2->two_point) ? "two-point" : "multipoint", 
	     (f1->two_point) ? "two-point" : "multipoint");
    exit (-1);
  }
  return;
}


/* compare_markers returns success or failure, instead of complaining and exiting like
 * the other comparison routines, because if a mapinfile has been provided, a marker
 * mismatch is not necessarily fatal.
 */
int compare_markers (st_brmarker *m1, st_brmarker *m2)
{
  if (! relax) {
    if (strcmp (m1->chr, m2->chr) != 0)
      return (-1);
    if ((strlen (m1->name1) != 0) && (strlen (m2->name1) != 0) &&
	(strcmp (m1->name1, m2->name1) != 0))
      return (-1);
    if ((strlen (m1->name2) != 0) && (strlen (m2->name2) != 0) &&
	(strcmp (m1->name2, m2->name2) != 0))
      return (-1);
  }
  return (0);
}


/* This behaves more like a 'strcmp' sort of affair, where it returns less than 0, 0,
 * or greater than 0 by comparing chromosome names and cM positions
 */
int compare_positions (st_brmarker *m1, st_brmarker *m2)
{
  int ret;
  double diff;

  if ((ret = strcmp (m1->chr, m2->chr)) != 0)
    return (ret);
  if ((diff = m1->pos - m2->pos) < 0)
    return (-1);
  else if (diff > 0)
    return (1);
  return (0);
}


void compare_samples (st_data *d, int sampleno, st_brfile *brfile)
{
  /* dcuhre2 imported from integrationLocals.h */
  if (fabs (dcuhre2[sampleno][0] - d->dprimes[0]) > 0.005) {
    fprintf (stderr, "D' varies from dKelvin at line %d in '%s'; expected %.2f, found %.2f\n",
	     brfile->lineno, brfile->name, dcuhre2[sampleno][0], d->dprimes[0]);
    exit (-1);
  }
  if (fabs (dcuhre2[sampleno][1] - d->thetas[0]) > 0.00005) {
    fprintf (stderr, "Theta varies from dKelvin at line %d in '%s'; expected %.4f, found %.4f\n",
	     brfile->lineno, brfile->name, dcuhre2[sampleno][1], d->thetas[0]);
    exit (-1);
  }
  
  /* I don't really think dKelvin supports sex-specific thetas... */
  if (! sexspecific) 
    return;
  if (fabs (dcuhre2[sampleno][1] - d->thetas[1]) > 0.00005) {
    fprintf (stderr, "Theta varies from dKelvin at line %d in '%s'; expected %.4f, found %.4f\n",
	     brfile->lineno, brfile->name, dcuhre2[sampleno][1], d->thetas[1]);
    exit (-1);
  }
  return;
}


#if 0
/* For sex-averaged thetas, we calculate the area under a two-dimensional
 * curve defined by theta (X axis) vs. LR (Y axis), where D' (all D' if
 * more than one) is 0. The area is calculated by computing the areas of
 * polygons bounded by adjascent thetas and associated LRs.
 */
double calc_ppl_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr)
{
  int zero_didx, thidx, va;
  double lr1, lr2, mtheta1, mtheta2, cutlr, pre_int=0, post_int=0, integral, ppl, *zeros;
  st_dim *mthetas;
  
  /* Remember for PPL, we only consider LRs where D' == 0, so look up that index
   * into the D' dimension of the lr array.
   */
  if ((zeros = malloc (sizeof (double) * dprimes->numdims)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes->numdims; va++)
    zeros[va] = 0.0;
  if ((zero_didx = multi_find (dprimes, zeros, dprimes->numdims)) == -1) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }
  free (zeros);

  /* if any of the LRs where D' == 0 is DBL_MAX (that is, undefined), then there
   * shouldn't be any defined LR where D' == 0, regardless of theta, and
   * we can't calculate a PPL.
   */
  if (lr[zero_didx][0] == DBL_MAX)
    return (-1);

  mthetas = &thetas->dims[0];
  for (thidx = 1; thidx < mthetas->numelems; thidx++) {
    mtheta1 = mthetas->arr[thidx-1];
    mtheta2 = mthetas->arr[thidx];
    lr1 = lr[zero_didx][thidx-1];
    lr2 = lr[zero_didx][thidx];

    if (mtheta1 >= cutoff) {
      /* entire region outside the cutoff */
      post_int += (lr1 + lr2) * (mtheta2 - mtheta1);
      
    } else if (mtheta2 <= cutoff) {
      /* entire region inside the cutoff */
      pre_int += (lr1 + lr2) * (mtheta2 - mtheta1);
      
    } else {
      /* region straddles the cutoff */
      cutlr = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
      pre_int += (lr1 + cutlr) * (cutoff - mtheta1);
      post_int += (cutlr + lr2) * (mtheta2 - cutoff);
    }
  }

  /* The 0.5 here is factored out of the calculation of area. */
  pre_int *= (weight / cutoff) * 0.5;
  post_int *= ((1 - weight) / (0.5 - cutoff)) * 0.5;
  integral = pre_int + post_int;
  ppl = (prior * integral) / (prior * integral + (1 - prior));
  return (ppl);
}


/* For sex-specific thetas, we calculate the volume under a three-dimensional
 * curve defined by male theta (X axis) vs. female theta (Y axis) vs. LR (Z
 * axis), where D' (all D', if more than one) is 0. The volume is calculated
 * by computing the volumes of polyhedrons bounded by adjascent males thetas,
 * adjascent female thetas, and the four LRs associated with the intersections
 * of the male and female thetas. Looking at the two-dimensional X-Y plane,
 * LR1 would be in the lower left, LR2 in the lower right, LR3 in the upper
 * left, and LR4 in the upper right.
 */
double calc_ppl_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr)
{
  int zero_didx, mthidx, fthidx, va;
  double lr1, lr2, lr3, lr4, mtheta1, mtheta2, ftheta1, ftheta2, cutlr2, cutlr3, cutlr4, *zeros;
  double pre_int=0, post_int=0, cutvol, integral, ppl;
  st_dim *mthetas, *fthetas;
  
  if ((zeros = malloc (sizeof (double) * dprimes->numdims)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprimes->numdims; va++)
    zeros[va] = 0.0;
  if ((zero_didx = multi_find (dprimes, zeros, dprimes->numdims)) == -1) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }
  free (zeros);

  /* See comment in calc_ppl_sexavg() */
  if (lr[zero_didx][0] == DBL_MAX)
    return (-1);
  
  /* Does it really matter which one is male or female? */
  mthetas = &thetas->dims[0];
  fthetas = &thetas->dims[1];
  for (mthidx = 1 ; mthidx < mthetas->numelems; mthidx++) {
    mtheta1 = mthetas->arr[mthidx - 1];
    mtheta2 = mthetas->arr[mthidx];
    for (fthidx = 1 ; fthidx < fthetas->numelems; fthidx++) {
      ftheta1 = fthetas->arr[fthidx - 1];
      ftheta2 = fthetas->arr[fthidx];
      lr1 = lr[zero_didx][(mthidx - 1) * mthetas->numelems + fthidx - 1];
      lr2 = lr[zero_didx][mthidx * mthetas->numelems + fthidx - 1];
      lr3 = lr[zero_didx][(mthidx - 1) * mthetas->numelems + fthidx];
      lr4 = lr[zero_didx][mthidx * mthetas->numelems + fthidx];
      
      if ((mtheta1 >= cutoff) || (ftheta1 >= cutoff)) {
	/* entire region is outside cutoff area */
	post_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	
      } else if ((mtheta2 <= cutoff) && (ftheta2 <= cutoff)) {
	/* entire region inside cutoff area */
	pre_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	
      } else if ((ftheta2 <= cutoff) && ((mtheta1 <= cutoff) && (cutoff <= mtheta2))) {
	/* region straddles cutoff on male axis */
	cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	pre_int += (cutoff - mtheta1) * (ftheta2 - ftheta1) * (lr1 + cutlr2 + lr3 + cutlr4);
	post_int += (mtheta2 - cutoff) * (ftheta2 - ftheta1) * (cutlr2 + lr2 + cutlr4 + lr4);
	
      } else if ((mtheta2 <= cutoff) && ((ftheta1 <= cutoff) && (cutoff <= ftheta2))) {
	/* region straddles cutoff on female axis */
	cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	cutlr4 = lr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr4 - lr2);
	pre_int += (cutoff - ftheta1) * (mtheta2 - mtheta1) * (lr1 + lr2 + cutlr3 + cutlr4);
	post_int += (ftheta2 - cutoff) * (mtheta2 - mtheta1) * (cutlr3 + cutlr4 + lr3 + lr4);
	
      } else {
	/* region straddles cutoff on both axes */
	post_int += (mtheta2 - mtheta1) * (ftheta2 - ftheta1) * (lr1 + lr2 + lr3 + lr4);
	cutlr2 = lr1 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr2 - lr1);
	cutlr4 = lr3 + ((cutoff - mtheta1) / (mtheta2 - mtheta1)) * (lr4 - lr3);
	cutlr4 = cutlr2 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (cutlr4 - cutlr2);
	cutlr3 = lr1 + ((cutoff - ftheta1) / (ftheta2 - ftheta1)) * (lr3 - lr1);
	cutvol = (cutoff - mtheta1) * (cutoff - ftheta1) * (lr1 + cutlr2 + cutlr3 + cutlr4);
	post_int -= cutvol;
	pre_int += cutvol;
      }
    }
  }

  /* The 0.25 here is factored out of the caclculation of volume */
  pre_int *= (weight / (cutoff * cutoff)) * 0.25;
  post_int *= ((1 - weight) / (0.25 - (cutoff * cutoff))) * 0.25;
  integral = pre_int + post_int;
  ppl = (prior * integral) / (prior * integral + (1 - prior));
  return (ppl);
}
#endif

