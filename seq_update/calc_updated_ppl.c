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

/* TO DO:
 * - sex-specific theta cutoffs
 * - multipoint PPLs
 * - mechanism for detecting missing D'/Theta combinations
 */

#define DPRIME_COL 1
#define THETA_COL  2
#define LR_COL     3
#define POS_COL    4
#define CHR_COL    5

#define METH_OLD   1
#define METH_NEW   2

#define BUFFLEN 256
#define KROUND(dbl) dbl >= 0.025 ? rint (dbl * 100.0) / 100.0 : rint (dbl * 10000.0) / 10000.0

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
    *datacols;
} st_brfile;

typedef struct {
  char name1[32],
    name2[32];
  double pos;
  int chr,
    num;
} st_marker;

typedef struct {
  int numdprimes,
    numthetas;
  double lr,
    *dprimes,
    *thetas;
} st_data;

typedef struct {
  double ld_small_theta,
    ld_big_theta,
    ld_unlinked,
    le_small_theta,
    le_big_theta,
    le_unlinked;
} st_ldvals;

typedef struct {
  int arrsize,
    numelems,
    lastidx,
    dimsize;
  double *arr;
} st_dim;

typedef struct {
  short numdims;
  int totalelems;
  st_dim *dims;
} st_multidim;


/* Default values that can be overridden on the command line */
int sexspecific = 0;    /* -s or --sexspecific: sex-specific thetas */
int multipoint = 0;     /* -m or --multipoint, not implemented yet */
int relax = 0;          /* -r or --relax: don't compare marker names between files */
double prior = 0.02;    /* -p or --prior: prior probability */
double ldprior = 0.021; /* -l or --ldprior: prior probability of LD */
double weight = 0.95;   /* -w or --weight: weighting factor */
double cutoff = 0.05;   /* -c or --cutoff: theta cutoff */
int method = METH_OLD;  /* --method: old (update BR per marker/theta/D' first) or 
			   new (integrate out theta/D' then update) */
FILE *partout = NULL;   /* --partout: file to which to write partial results */
FILE *partin = NULL;    /* --partin: file from which to read partial ldval-style fields */

/* Globally available for error messages */
char *current = "0.36.1";
char *pname;
char *partinfile=NULL;
/* Globbaly available for building partial output files */
char partheader[BUFFLEN], partpadding[BUFFLEN];

void old_method (st_brfile *brfiles, int numbrfiles);
void old_method_multipoint (st_brfile *brfiles, int numbrfiles);
void new_method (st_brfile *brfiles, int numbrfiles);
void do_first_pass (st_brfile *brfile, st_multidim *dprimes, st_multidim *thetas, st_data *data);
double calc_ppl_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr);
double calc_ppl_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr);
void calc_ldvals_sexavg (st_multidim *dprimes, st_multidim *thetas, double **lr, st_ldvals *ldvals);
void calc_ldvals_sexspc (st_multidim *dprimes, st_multidim *thetas, double **lr, st_ldvals *ldvals);
double calc_ldppl (st_ldvals *ldvals);
double calc_ppld_given_linkage (st_ldvals *ldvals);
double calc_ppld (st_ldvals *ldvals);
double calc_ppld_and_linkage (st_ldvals *ldvals);

int parse_command_line (int argc, char **argv);
void usage ();
double validate_double_arg (char *arg, char *optname);
void open_brfile (st_brfile *brfile);
int get_marker_line (st_brfile *brfile, st_marker *marker);
int get_header_line (st_brfile *brfile, st_data *data);
int get_data_line (st_brfile *brfile, st_marker *marker, st_data *data);
void read_partin (st_marker *markers, st_ldvals *ldvals, int nummarkers);
void print_partial_old (st_brfile *brfile, st_marker *marker, st_data *data);
void compare_headers (st_brfile *f1, st_brfile *f2);
void compare_markers (st_marker *m1, st_marker *m2, st_brfile *brfile);
void compare_positions (st_marker *m1, st_marker *m2, st_brfile *brfile);
int multi_insert (st_multidim *md, double *vals, int num);
int multi_find (st_multidim *md, double *vals, int num);
int insert (st_dim *dim, double val);
int find (st_dim *dim, double val, int *found);
void multi_free (st_multidim *md);
void multi_dump (st_multidim *md);


int main (int argc, char **argv)
{
  int argidx, va, numbrfiles;
  st_brfile *brfiles;

  pname = argv[0];
  partheader[0] = partpadding[0] = '\0';
  argidx = parse_command_line (argc, argv);
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

  if (method == METH_OLD) {
    if (multipoint) {
      old_method_multipoint (brfiles, numbrfiles);
    } else {
      old_method (brfiles, numbrfiles);
    }
  } else {
    new_method (brfiles, numbrfiles);
  }
  for (va = 0; va < numbrfiles; va++) {
    fclose (brfiles[va].fp);
  }
  free (brfiles);
  if (partout != NULL)
    fclose (partout);
  exit (0);
}


void old_method (st_brfile *brfiles, int numbrfiles)
{
  int ret, fileno, didx, thidx, firstpass = 1;
  double **lr, ldstat;
  st_marker marker_0, marker_n;
  st_data data;
  st_multidim dprimes, thetas;
  st_ldvals ldval;

  memset (&marker_0, 0, sizeof (st_marker));
  memset (&marker_n, 0, sizeof (st_marker));
  memset (&data, 0, sizeof (st_data));
  memset (&dprimes, 0, sizeof (st_multidim));
  memset (&thetas, 0, sizeof (st_multidim));

  /* First thing to do is read through the set of data lines for the first
   * marker in the first file. We can use the D' and thetas that we read to 
   * populate our data structures. Since this is the old method, every file
   * has to have the same D' and thetas for every marker.
   */
  do_first_pass (&brfiles[0], &dprimes, &thetas, &data);

  if ((lr = malloc (sizeof (double *) * dprimes.totalelems)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (didx = 0; didx < dprimes.totalelems; didx++) {
    if ((lr[didx] = malloc (sizeof (double) * thetas.totalelems)) == NULL) {
      fprintf (stderr, "malloc failed, %s\n", strerror (errno));
      exit (-1);
    }
  }

  if (! brfiles[0].no_ld)
    printf ("Chr Seq Marker Position PPL LD-PPL PPLD|L PPLD PPLD&L\n");
  else
    printf ("Chr Seq Marker Position PPL\n");


  while ((ret = get_marker_line (&brfiles[0], &marker_0)) == 1) {
    get_header_line (&brfiles[0], &data);
    for (fileno = 1; fileno < numbrfiles; fileno++) {
      if (get_marker_line (&brfiles[fileno], &marker_n) != 1) {
	fprintf (stderr, "file '%s' ends where marker line expected at line %d\n",
		 brfiles[fileno].name, brfiles[fileno].lineno);
	exit (-1);
      }
      get_header_line (&brfiles[fileno], &data);
      if (firstpass)
	compare_headers (&brfiles[0], &brfiles[fileno]);
      compare_markers (&marker_0, &marker_n, &brfiles[fileno]);
    }   
    firstpass = 0;
    if (partout != NULL) 
      fprintf (partout, "# %d %s %s\n%s", marker_0.num, marker_0.name1, marker_0.name2, 
	       partheader);
    
    for (didx = 0; didx < dprimes.totalelems; didx++) {
      for (thidx = 0; thidx < thetas.totalelems; thidx++) {
	lr[didx][thidx] = DBL_MAX;
      }
    }

    while ((ret = get_data_line (&brfiles[0], &marker_0, &data)) == 1) {
      if ((didx = multi_find (&dprimes, data.dprimes, brfiles[0].numdprimes)) == -1) {
	fprintf (stderr, "unexpected dprimes in '%s' at line %d\n", brfiles[0].name,
		 brfiles[0].lineno);
	exit (-1);
      }
      if ((thidx = multi_find (&thetas, data.thetas, brfiles[0].numthetas)) == -1) {
	fprintf (stderr, "unexpected thetas in '%s' at line %d\n", brfiles[0].name,
		 brfiles[0].lineno);
	exit (-1);
      }
      if (lr[didx][thidx] != DBL_MAX) {
	fprintf (stderr, "duplicate dprime/theta combination in '%s' at line %d\n",
		 brfiles[0].name, brfiles[0].lineno);
	exit (-1);
      }
      lr[didx][thidx] = data.lr;
      
      for (fileno = 1; fileno < numbrfiles; fileno++) {
	if ((ret = get_data_line (&brfiles[fileno], &marker_n, &data)) == -1) {
	  fprintf (stderr, "can't parse data, line %d in file '%s'\n", brfiles[fileno].lineno,
		   brfiles[fileno].name);
	} else if (ret != 1) {
	  fprintf (stderr, "data ends unexpectedly at line %d in file '%s'\n",
		   brfiles[fileno].lineno, brfiles[fileno].name);
	}
	if (multi_find (&dprimes, data.dprimes, brfiles[fileno].numdprimes) != didx) {
	  fprintf (stderr, "unexpected dprimes in '%s' at line %d\n", brfiles[fileno].name,
		   brfiles[fileno].lineno);
	  exit (-1);
	}
	if (multi_find (&thetas, data.thetas, brfiles[fileno].numthetas) != thidx) {
	  fprintf (stderr, "unexpected thetas in '%s' at line %d\n", brfiles[fileno].name,
		   brfiles[fileno].lineno);
	  exit (-1);
	}
	lr[didx][thidx] *= data.lr;
      }
      data.lr = lr[didx][thidx];
      if (partout != NULL)
	print_partial_old (&brfiles[0], &marker_0, &data);
    }
    if (ret == -1) {
      fprintf (stderr, "can't parse data, line %d in file '%s'\n", brfiles[0].lineno,
	       brfiles[0].name);
      exit (-1);
    }

    printf ("%d %d %s %.4f", marker_0.chr, marker_0.num, marker_0.name2, marker_0.pos);
    printf (" %.3f", (! sexspecific) ? calc_ppl_sexavg (&dprimes, &thetas, lr) : 
	    calc_ppl_sexspc (&dprimes, &thetas, lr));
    if (! brfiles[0].no_ld) {
      if (! sexspecific)
	calc_ldvals_sexavg (&dprimes, &thetas, lr, &ldval);
      else 
	calc_ldvals_sexspc (&dprimes, &thetas, lr, &ldval);
      
      ldstat = calc_ldppl (&ldval);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_ppld_given_linkage (&ldval);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_ppld (&ldval);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_ppld_and_linkage (&ldval);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
    }
    printf ("\n");
  }
  if (ret == -1) {
    fprintf (stderr, "can't parse marker, line %d in file '%s'\n", brfiles[0].lineno,
	     brfiles[0].name);
    exit (-1);
  }
  
  for (didx = 0; didx < dprimes.totalelems; didx++)
    free (lr[didx]);
  free (lr);

  multi_free (&dprimes);
  multi_free (&thetas);
  
  free (data.dprimes);
  free (data.thetas);
  for (fileno = 0; fileno < numbrfiles; fileno++)
    free (brfiles[fileno].datacols);
  return;
}


void old_method_multipoint (st_brfile *brfiles, int numbrfiles)
{
  int fileno;
  st_marker marker_0, marker_n;
  st_data data;
  double lr, ppl;

  memset (&data, 0, sizeof (st_data));

  get_header_line (&brfiles[0], &data);
  for (fileno = 1; fileno < numbrfiles; fileno++) {
    get_header_line (&brfiles[fileno], &data);
    compare_headers (&brfiles[0], &brfiles[fileno]);
  }
  printf ("Chr Position PPL BayesRatio\n");

  while (get_data_line (&brfiles[0], &marker_0, &data) == 1) {
    lr = data.lr;
    for (fileno = 1; fileno < numbrfiles; fileno++) {
      if (get_data_line (&brfiles[fileno], &marker_n, &data) != 1) {
	fprintf (stderr, "file '%s' ends unexpectedly at line %d\n",
		 brfiles[fileno].name, brfiles[fileno].lineno);
	exit (-1);
      }
       compare_positions (&marker_0, &marker_n, &brfiles[fileno]);
      lr *= data.lr;
    }
    if ((lr < 0.214) || ((ppl = (lr * lr) / (-5.77 + (54 * lr) + (lr * lr))) < 0.0))
      ppl = 0.0;
    printf ("%d %.4f %.2f %.6e\n", marker_0.chr, marker_0.pos, ppl, lr);
  }
  return;
}


void new_method (st_brfile *brfiles, int numbrfiles)
{
  int ret, fileno, mrkno, nummarkers=0, didx, thidx;
  double **lr, integral, ldstat;
  st_marker *markers=NULL, marker;
  st_ldvals *ldvals=NULL, ldval;
  st_multidim dprimes, thetas;
  st_data data;

  memset (&data, 0, sizeof (st_data));
  for (fileno = 0; fileno < numbrfiles; fileno++) {
    memset (&marker, 0, sizeof (st_marker));
    memset (&dprimes, 0, sizeof (st_multidim));
    memset (&thetas, 0, sizeof (st_multidim));
    mrkno = -1;
    
    do_first_pass (&brfiles[fileno], &dprimes, &thetas, &data);

    if ((lr = malloc (sizeof (double *) * dprimes.totalelems)) == NULL) {
      fprintf (stderr, "malloc failed, %s\n", strerror (errno));
      exit (-1);
    }
    for (didx = 0; didx < dprimes.totalelems; didx++) {
      if ((lr[didx] = malloc (sizeof (double) * thetas.totalelems)) == NULL) {
	fprintf (stderr, "malloc failed, %s\n", strerror (errno));
	exit (-1);
      }
    }

    while ((ret = get_marker_line (&brfiles[fileno], &marker)) == 1) {
      if (++mrkno >= nummarkers) {
	if (fileno == 0) {
	  if (((markers = realloc (markers, sizeof (st_marker) * ++nummarkers)) == NULL) ||
	      ((ldvals = realloc (ldvals, sizeof (st_ldvals) * nummarkers)) == NULL)) {
	    fprintf (stderr, "realloc failed, %s\n", strerror (errno));
	    exit (-1);
	  }
	} else {
	  fprintf (stderr, "file '%s' contains data for more markers than '%s'\n", 
		   brfiles[fileno].name, brfiles[0].name);
	  exit (-1);
	}
      }
      get_header_line (&brfiles[fileno], &data);
      if ((fileno != 0) && (brfiles[0].no_ld != brfiles[fileno].no_ld)) {
	fprintf (stderr, "can't mix %s ('%s') and %s ('%s') analyses\n",
		 (brfiles[0].no_ld) ? "non-LD" : "LD", brfiles[0].name,
		 (brfiles[fileno].no_ld) ? "non-LD" : "LD", brfiles[fileno].name);
	exit (-1);
      }
      
      for (didx = 0; didx < dprimes.totalelems; didx++) {
	for (thidx = 0; thidx < thetas.totalelems; thidx++) {
	  lr[didx][thidx] = DBL_MAX;
	}
      }
      
      while ((ret = get_data_line (&brfiles[fileno], &marker, &data)) == 1) {
	if ((didx = multi_find (&dprimes, data.dprimes, brfiles[fileno].numdprimes)) == -1) {
	  fprintf (stderr, "unexpected dprimes in '%s' at line %d\n", brfiles[fileno].name,
		   brfiles[fileno].lineno);
	  exit (-1);
	}
	if ((thidx = multi_find (&thetas, data.thetas, brfiles[fileno].numthetas)) == -1) {
	  fprintf (stderr, "unexpected thetas in '%s' at line %d\n", brfiles[fileno].name,
		   brfiles[fileno].lineno);
	  exit (-1);
	}
	if (lr[didx][thidx] != DBL_MAX) {
	  fprintf (stderr, "duplicate dprime/theta combination in '%s' at line %d\n",
		   brfiles[fileno].name, brfiles[fileno].lineno);
	  exit (-1);
	}
	lr[didx][thidx] = data.lr;
      }
      if (! sexspecific)
	calc_ldvals_sexavg (&dprimes, &thetas, lr, &ldval);
      else 
	calc_ldvals_sexspc (&dprimes, &thetas, lr, &ldval);
      
      if (fileno == 0) {
	memcpy (&markers[mrkno], &marker, sizeof (st_marker));
	memcpy (&ldvals[mrkno], &ldval, sizeof (st_ldvals));
      } else {
	compare_markers (&markers[mrkno], &marker, &brfiles[fileno]);
	ldvals[mrkno].ld_small_theta *= ldval.ld_small_theta;
	ldvals[mrkno].ld_big_theta *= ldval.ld_big_theta;
	ldvals[mrkno].ld_unlinked *= ldval.ld_unlinked;
	ldvals[mrkno].le_small_theta *= ldval.le_small_theta;
	ldvals[mrkno].le_big_theta *= ldval.le_big_theta;
	ldvals[mrkno].le_unlinked *= ldval.le_unlinked;
      }
    }
    
    for (didx = 0; didx < dprimes.totalelems; didx++) {
      free (lr[didx]);
    }
    free (lr);
    
    multi_free (&dprimes);
    multi_free (&thetas);
    
    free (brfiles[fileno].datacols);
  }

  if (partin != NULL)
    read_partin (markers, ldvals, nummarkers);

  if (! brfiles[0].no_ld)
    printf ("Chr Seq Marker Position PPL LD-PPL PPLD|L PPLD PPLD&L\n");
  else
    printf ("Chr Seq Marker Position PPL\n");

  if (partout != NULL)
    fprintf (partout, "Chr Pos Seq Name1 Name2 LDSmallTheta LDBigTheta LDUnlinked LESmallTheta LEBigTheta LEUnlinked\n");
  
  for (mrkno = 0; mrkno < nummarkers; mrkno++) {
    printf ("%d %d %s %.4f", markers[mrkno].chr, markers[mrkno].num, markers[mrkno].name2,
	    markers[mrkno].pos);

    /* print PPL - is this right? At least it takes up space... */
    integral = ldvals[mrkno].le_small_theta + ldvals[mrkno].le_big_theta;
    printf (" %.3f", (prior * integral) / (prior * integral + (1 - prior)));
    
    if (! brfiles[0].no_ld) {
      ldstat = calc_ldppl (&ldvals[mrkno]);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_ppld_given_linkage (&ldvals[mrkno]);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_ppld (&ldvals[mrkno]);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
      ldstat = calc_ppld_and_linkage (&ldvals[mrkno]);
      printf (" %.*f", ldstat >= .025 ? 2 : 4, KROUND (ldstat));
    }
    printf ("\n");

    if (partout != NULL)
      fprintf (partout, "%d %.4f %d %s %s %.6e %.6e %.6e %.6e %.6e %.6e\n", markers[mrkno].chr,
	       markers[mrkno].pos, markers[mrkno].num, markers[mrkno].name1, markers[mrkno].name2,
	       ldvals[mrkno].ld_small_theta, ldvals[mrkno].ld_big_theta,
	       ldvals[mrkno].ld_unlinked, ldvals[mrkno].le_small_theta,
	       ldvals[mrkno].le_big_theta, ldvals[mrkno].le_unlinked);

  }
  free (data.dprimes);
  free (data.thetas);
  free (markers);
  free (ldvals);
  return;
}


void do_first_pass (st_brfile *brfile, st_multidim *dprimes, st_multidim *thetas, st_data *data)
{
  int ret;
  long firstmarker;
  st_marker marker;

  memset (&marker, 0, sizeof (st_marker));

  if ((firstmarker = ftell (brfile->fp)) == -1) {
    fprintf (stderr, "ftell on file '%s' failed, %s\n", brfile->name, strerror (errno));
    exit (-1);
  }
  if (get_marker_line (brfile, &marker) == 0) {
    fprintf (stderr, "file '%s' is empty\n", brfile->name);
    exit (-1);
  }
  get_header_line (brfile, data);

  while ((ret = get_data_line (brfile, &marker, data)) == 1) {
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
  
  if (fseek (brfile->fp, firstmarker, SEEK_SET) == -1) {
    fprintf (stderr, "fseek on file '%s' failed, %s\n", brfile->name, strerror (errno));
    exit (-1);
  }
  brfile->lineno = 1;
  return;
}


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


double calc_ldppl (st_ldvals *ldval)
{
  double numerator;
  double denomRight;
  double ldppl;

  numerator =
    ldval->ld_small_theta * prior * 0.021 + 
    ldval->ld_big_theta * prior *  0.0011+ 
    ldval->le_small_theta * prior * 0.979 + 
    ldval->le_big_theta * prior * 0.9989;
  denomRight =
    ldval->le_unlinked * (1 - prior);
  ldppl = numerator / (numerator + denomRight);
  
  return (ldppl);
}


double calc_ppld (st_ldvals *ldval)
{
  double numerator;
  double denomRight;
  double ppld; 

  numerator =
    ldval->ld_small_theta * prior * 0.021 +
    ldval->ld_big_theta * prior * 0.0011; 
  denomRight =
    ldval->le_small_theta * prior * 0.979 + 
    ldval->le_big_theta * prior * 0.9989 + 
    ldval->le_unlinked * (1 - prior);
  ppld = numerator / (numerator + denomRight); 

  return (ppld);
}


double calc_ppld_and_linkage (st_ldvals *ldval)
{
  double numerator;
  double denomRight;
  double ppld_and_l;

  numerator =
    ldval->ld_small_theta * prior * 0.021 + 
    ldval->ld_big_theta * prior * 0.0011;
  denomRight =
    ldval->le_small_theta * prior * 0.979 + 
    ldval->le_big_theta * prior * 0.9989 + 
    ldval->le_unlinked * (1 - prior);
  ppld_and_l = numerator / (numerator + denomRight);

  return (ppld_and_l);
}


double calc_ppld_given_linkage (st_ldvals *ldval)
{
  double numerator;
  double denomRight;
  double ppld_given_l;

  numerator =
    ldval->ld_small_theta * prior * 0.021 + 
    ldval->ld_big_theta * prior * 0.0011;
  denomRight =
    ldval->le_small_theta * prior * 0.979 + 
    ldval->le_big_theta * prior * 0.9989;
  ppld_given_l = numerator / (numerator + denomRight);

  return (ppld_given_l);
}


#define OPT_SEXSPEC 1
#define OPT_MULTI   2
#define OPT_RELAX   3
#define OPT_PRIOR   4
#define OPT_LDPRIOR 5
#define OPT_WEIGHT  6
#define OPT_CUTOFF  7
#define OPT_METHOD  8
#define OPT_PARTIN  9
#define OPT_PARTOUT 10
#define OPT_HELP    11

int parse_command_line (int argc, char **argv)
{
  int arg, long_arg, long_idx;
  char *partoutfile=NULL;
  struct option cmdline[] = { { "sexspecific", 0, &long_arg, OPT_SEXSPEC },
			      { "multipoint", 0, &long_arg, OPT_MULTI },
			      { "relax", 1, &long_arg, OPT_RELAX },
			      { "prior", 1, &long_arg, OPT_PRIOR },
			      { "ldprior", 1, &long_arg, OPT_LDPRIOR},
			      { "weight", 1, &long_arg, OPT_WEIGHT },
			      { "cutoff", 1, &long_arg, OPT_CUTOFF },
			      { "method", 1, &long_arg, OPT_METHOD },
			      { "partin", 1, &long_arg, OPT_PARTIN },
			      { "partout", 1, &long_arg, OPT_PARTOUT },
			      { "help", 0, &long_arg, OPT_HELP },
			      { NULL, 0, NULL, 0 } };
  struct stat statbuf;
  
  while ((arg = getopt_long (argc, argv, "smrp:l:w:c:", cmdline, &long_idx)) != -1) {
    if ((arg == 's') || ((arg == 0) && (long_arg == OPT_SEXSPEC))) {
      sexspecific = 1;

    } else if ((arg == 'm') || ((arg == 0) && (long_arg == OPT_MULTI))) {
      multipoint = 1;

    } else if ((arg == 'r') || ((arg == 0) && (long_arg == OPT_RELAX))) {
      relax = 1;

    } else if (arg == 'p') {
      prior = validate_double_arg (optarg, "-p");
    } else if ((arg == 0) && (long_arg == OPT_PRIOR)) {
      prior = validate_double_arg (optarg, "--prior");

    } else if (arg == 'l') {
      ldprior = validate_double_arg (optarg, "-l");
    } else if ((arg == 0) && (long_arg == OPT_LDPRIOR)) {
      ldprior = validate_double_arg (optarg, "--ldprior");

    } else if (arg == 'w') {
      weight = validate_double_arg (optarg, "-w");
    } else if ((arg == 0) && (long_arg == OPT_WEIGHT)) {
      weight = validate_double_arg (optarg, "--weight");

    } else if (arg == 'c') {
      cutoff = validate_double_arg (optarg, "-c");
    } else if ((arg == 0) && (long_arg == OPT_CUTOFF)) {
      cutoff = validate_double_arg (optarg, "--cutoff");

    } else if ((arg == 0) && (long_arg == OPT_METHOD)) {
      if (strcasecmp (optarg, "new") == 0) {
	method = METH_NEW;
      } else if (strcasecmp (optarg, "old") != 0) {
	fprintf (stderr, "%s: argument '%s' to option '--method' is illegal\n", pname, optarg);
	exit (-1);
      }

    } else if ((arg == 0) && (long_arg == OPT_PARTIN)) {
      partinfile = optarg;

    } else if ((arg == 0) && (long_arg == OPT_PARTOUT)) {
      partoutfile = optarg;

    } else if ((arg == 0) && (long_arg == OPT_HELP)) {
      usage ();
      exit (0);

    } else {
      /* getopt_long will emit a message describing the bad option/argument */
      usage ();
      exit (-1);
    }
  }

  if ((partinfile != NULL) && (method != METH_NEW)) {
    fprintf (stderr, "%s: don't use --partin without --method=new\n", pname);
    exit (-1);
  }

  if ((multipoint) && (method == METH_NEW)) {
    fprintf (stderr, "%s: --multipoint is nonsensical with --method=new\n", pname);
    exit (-1);
  }

  if (partinfile != NULL) {
    if ((partin = fopen (partinfile, "r")) == NULL) {
      fprintf (stderr, "%s: open '%s' for reading failed, %s\n", pname, partinfile,
	       strerror (errno));
      exit (-1);
    }
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
  return (optind);
}


void usage ()
{
  printf ("usage: %s [ options ] brfile [brfile...]\n  valid options:\n", pname);
  printf ("  -s|--sexspecific : input data contains sex-specific Thetas\n");
  printf ("  -m|--multipoint : input data is multipoint\n");
  printf ("  -r|--relax : supress comparing marker names across input files\n");
  printf ("  -p <num>|--prior <num> : set linkage prior probability to <num>\n");
  printf ("  -c <num>|--cutoff <num> : set small-Theta cutoff to <num>\n");
  printf ("  -w <num>|--wieght <num> : set small-Theta weight to <num>\n");
  printf ("  --partout <file> : write updated Bayes Ratios to <file>\n");
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
    if (feof (brfile->fp))
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


int get_marker_line (st_brfile *brfile, st_marker *marker)
{
  char buff[BUFFLEN];
  char *pa, *pb;
  
  if (fgets (buff, BUFFLEN, brfile->fp) == NULL) {
    /* An end-of-file looking for a marker is not necessarily fatal, so let the caller decide */
    if (feof (brfile->fp))
      return (0);
    fprintf (stderr, "error reading '%s' at line %d, %s\n", brfile->name, brfile->lineno,
	     strerror (errno));
    exit (-1);
  }
  brfile->lineno++;

  if (((pa = strtok_r (buff, " \t\n", &pb)) == NULL) ||
      (strcmp (pa, "#") != 0)) {
    fprintf (stderr, "can't parse marker, line %d in file '%s'\n", brfile->lineno,  brfile->name);
    exit (-1);
  }
  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL) {
    fprintf (stderr, "can't parse marker, line %d in file '%s'\n", brfile->lineno,  brfile->name);
    exit (-1);
  }
  if (((marker->num = (int) strtol (pa, NULL, 10)) == 0) && (errno != 0)) {
    fprintf (stderr, "can't parse marker, line %d in file '%s'\n", brfile->lineno,  brfile->name);
    exit (-1);
  }
  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL) {
    fprintf (stderr, "can't parse marker, line %d in file '%s'\n", brfile->lineno,  brfile->name);
    exit (-1);
  }
  strcpy (marker->name1, pa);
  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL) {
    fprintf (stderr, "can't parse marker, line %d in file '%s'\n", brfile->lineno,  brfile->name);
    exit (-1);
  }
  strcpy (marker->name2, pa);
  
  return (1);
}


int get_header_line (st_brfile *brfile, st_data *data)
{ 
  char buff[BUFFLEN], token[32], tmppadding[BUFFLEN];
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

  /* If datacols is non-NULL (that is, memory has been allocated) then we've parsed
   * a header for this file before. Since we require files to be internally consistent,
   * we can just return.
   */
  if (brfile->datacols != NULL)
    return (1);
  if ((method == METH_OLD) && (partout != NULL) && (strlen (partheader) == 0)) {
    strcpy (partheader, buff);
    tmppadding[0] = '\0';
  }
      
  pa = strtok_r (buff, " \t\n", &pb);
  while (pa != NULL) {
#ifdef DEBUG
    printf ("initial token '%s'\n", pa);
#endif
    strcpy (token, pa);
    while ((strchr (token, '(') != NULL) && (strrchr (token, ')') == NULL)) {
      if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL) {
	fprintf (stderr, "can't parse header, line %d in file '%s'\n", brfile->lineno,
		 brfile->name);
	exit (-1);
      }
#ifdef DEBUG
      printf ("  tacking on '%s'\n", pa);
#endif
      strcat (token, pa);
    }
#ifdef DEBUG
    printf ("final token '%s'\n", token);
#endif
    
    if (token[0] == '(') {
      fprintf (stderr, "i don't know why this conditional is here\n");
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
#ifdef DEBUG
    printf ("token is '%s', actualcols %d\n", token, actualcols);
#endif

    brfile->numcols += actualcols;
#ifdef DEBUG
    printf ("reallocating datacols to %d\n", brfile->numcols);
#endif
    if ((brfile->datacols = realloc (brfile->datacols, sizeof (int) * brfile->numcols)) == NULL) {
      fprintf (stderr, "realloc failed, %s\n", strerror (errno));
      exit (-1);
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
      if ((method == METH_OLD) && (partout != NULL) && (strlen (partpadding) == 0)) {
	if (actualcols == 1) {
	  strcat (tmppadding, " 0");
	} else {
	  strcat (tmppadding, " (0");
	  for (va = actualcols-1; va > 0; va--) 
	    strcat (tmppadding, ",0");
	  strcat (tmppadding, ")");
	}
      }
      for (; actualcols > 0; actualcols--) 
        brfile->datacols[brfile->numcols - actualcols] = 0;
    }

#ifdef DEBUG
    printf ("\n");
#endif
    pa = strtok_r (NULL, " \t\n", &pb);
  }

  if ((method == METH_OLD) && (partout != NULL) && (strlen (partpadding) == 0))
    strcpy (partpadding, tmppadding);

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
  
  if ((data->numdprimes < brfile->numdprimes) && 
      ((data->dprimes = realloc (data->dprimes, sizeof (double) * brfile->numdprimes)) == NULL)) {
    fprintf (stderr, "realloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  data->numdprimes = brfile->numdprimes;
  if ((data->numthetas < brfile->numthetas) && 
      ((data->thetas = realloc (data->thetas, sizeof (double) * brfile->numthetas)) == NULL)) {
    fprintf (stderr, "realloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  data->numthetas = brfile->numthetas;

  #ifdef DEBUG
  for (numlrcols = 0; numlrcols < brfile->numcols; numlrcols++) {
    printf ("col %d is type %d\n", numlrcols, brfile->datacols[numlrcols]);
  }
  #endif

  return (1);
}


int get_data_line (st_brfile *brfile, st_marker *marker, st_data *data)
{
  char buff[BUFFLEN], *endptr=NULL;
  char *pa=NULL, *pb=NULL;
  long previous;
  int va = 0, dprimecnt=0, thetacnt=0;

  if ((previous = ftell (brfile->fp)) == -1) {
    fprintf (stderr, "ftell failed, %s\n", strerror (errno));
    exit (-1);
  }
  if (fgets (buff, BUFFLEN, brfile->fp) == NULL) {
    if (feof (brfile->fp))
      return (0);
    fprintf (stderr, "error reading '%s' at line %d, %s\n", brfile->name, brfile->lineno,
	     strerror (errno));
    exit (-1);
  }
  brfile->lineno++;
  
  if ((pa = strtok_r (buff, " (),\t\n", &pb)) == NULL) {
    fprintf (stderr, "unexpectedly short line %d in '%s'\n", brfile->lineno, brfile->name);
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
      marker->chr = (int) strtol (pa, &endptr, 10);
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
  printf ("%5d: chr %2d, pos %6.4f, dprimes", brfile->lineno, marker->chr, marker->pos);
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


void read_partin (st_marker *markers, st_ldvals *ldvals, int nummarkers)
{
  int mrkno=0;
  char buff[BUFFLEN];
  st_brfile partbr;
  st_marker marker;
  st_ldvals ldval;

  if (fgets (buff, BUFFLEN, partin) == NULL) {
    if (feof (partin))
      fprintf (stderr, "partin file '%s' is empty\n", partinfile);
    else
      fprintf (stderr, "error reading partin file '%s', %s\n", partinfile, strerror (errno));
    exit (-1);
  }
  /* skip detailed version verification for the time being */
  if (strstr (buff, "Version") == NULL) {
    fprintf (stderr, "file %s is unversioned\n", partinfile);
    exit (-1);
  }
  
  if (fgets (buff, BUFFLEN, partin) == NULL) {
    if (feof (partin))
      fprintf (stderr, "partin file '%s' is ends unexpectedly at line 1\n", partinfile);
    else
      fprintf (stderr, "error reading partin file '%s', %s\n", partinfile, strerror (errno));
    exit (-1);
  }
  if (strstr (buff, "Chr Pos Seq Name1 Name2 LDSmallTheta") == NULL) {
    fprintf (stderr, "bad header line in '%s' at line 2\n", partinfile);
    exit (-1);
  }
  
  partbr.name = partinfile;
  partbr.lineno = 2;

  for (mrkno = 0; mrkno < nummarkers; mrkno++) {
    if (fgets (buff, BUFFLEN, partin) == NULL) {
      if (feof (partin))
	fprintf (stderr, "partin file '%s' is ends unexpectedly at line 1\n", partinfile);
      else
	fprintf (stderr, "error reading partin file '%s', %s\n", partinfile, strerror (errno));
      exit (-1);
    }
    partbr.lineno++;
    if (sscanf (buff, "%d %lf %d %s %s %lf %lf %lf %lf %lf %lf",
		&marker.chr, &marker.pos, &marker.num, marker.name1, marker.name2,
		&ldval.ld_small_theta, &ldval.ld_big_theta, &ldval.ld_unlinked,
		&ldval.le_small_theta, &ldval.le_big_theta, &ldval.le_unlinked) != 11) {
      fprintf (stderr, "can't parse line %d in '%s'\n", partbr.lineno, partinfile);
      exit (-1);
    }
    compare_markers (&markers[mrkno], &marker, &partbr);
    ldvals[mrkno].ld_small_theta *= ldval.ld_small_theta;
    ldvals[mrkno].ld_big_theta *= ldval.ld_big_theta;
    ldvals[mrkno].ld_unlinked *= ldval.ld_unlinked;
    ldvals[mrkno].le_small_theta *= ldval.le_small_theta;
    ldvals[mrkno].le_big_theta *= ldval.le_big_theta;
    ldvals[mrkno].le_unlinked *= ldval.le_unlinked;
  }
  fclose (partin);
  return;
}


void print_partial_old (st_brfile *brfile, st_marker *marker, st_data *data)
{
  int colno=0, dprimecnt=0;

  while (colno < brfile->numcols) {
    switch (brfile->datacols[colno]) {
    case CHR_COL:
      fprintf (partout, "%d", marker->chr);
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
  fprintf (partout, "%s\n", partpadding);
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


void compare_markers (st_marker *m1, st_marker *m2, st_brfile *brfile)
{
  if (! relax) {
    if (m1->num != m2->num) {
      fprintf (stderr, "marker number mismatch at line %d in '%s'; expected %d, found %d\n",
	       brfile->lineno, brfile->name, m1->num, m2->num);
      exit (-1);
    }
    if (strcmp (m1->name1, m2->name1) != 0) {
      fprintf (stderr, "marker name mismatch at line %d in '%s'; expected %s, found %s\n",
	       brfile->lineno, brfile->name, m1->name1, m2->name1);
      exit (-1);
    }
    if (strcmp (m1->name2, m2->name2) != 0) {
      fprintf (stderr, "marker name mismatch at line %d in '%s'; expected %s, found %s\n",
	       brfile->lineno, brfile->name, m1->name2, m2->name2);
      exit (-1);
    }
  }
  return;
}


void compare_positions (st_marker *m1, st_marker *m2, st_brfile *brfile)
{
  if (! relax) {
    if (m1->chr != m2->chr) {
      fprintf (stderr, "chromosome number mismatch at line %d in '%s'; expected %d, found %d\n",
	       brfile->lineno, brfile->name, m1->chr, m2->chr);
      exit (-1);
    }
    if (m1->pos != m2->pos) {
      fprintf (stderr, "position mismatch at line %d in '%s'; expected %f, found %f\n",
	       brfile->lineno, brfile->name, m1->pos, m2->pos);
      exit (-1);
    }
  }
  return;
}




int multi_insert (st_multidim *md, double *vals, int num)
{
  int va, idx, offset=0;
  st_dim *tmp;

  if (md->numdims != num) {
    if ((tmp = realloc (md->dims, sizeof (st_dim) * num)) == NULL)
      return (-1);
    md->dims = tmp;
    for (va = md->numdims; va < num; va++)
      memset (&(md->dims[va]), 0, sizeof (st_dim));
    md->numdims = num;
  }
  
  for (va = md->numdims - 1; va >= 0; va--) {
    if ((idx = insert (&(md->dims[va]), vals[va])) == -1)
      return (-1);
    if (va == md->numdims - 1)
      md->dims[va].dimsize = 1;
    else
      md->dims[va].dimsize = md->dims[va+1].numelems * md->dims[va+1].dimsize;
    offset += md->dims[va].dimsize * idx;
  }
  md->totalelems = md->dims[0].dimsize * md->dims[0].numelems;
  return (offset);
}


int multi_find (st_multidim *md, double *vals, int num)
{
  int va, idx, found, offset=0;

  if (md->numdims != num)
    return (-1);
  for (va = md->numdims - 1; va >= 0; va--) {
    idx = find (&(md->dims[va]), vals[va], &found);
    if (! found)
      return (-1);
    offset += idx * md->dims[va].dimsize;
  }
  return (offset);
}


int insert (st_dim *dim, double val)
{
  double *tmp;
  int idx, found;

  idx = find (dim, val, &found);
  if (found)
    return (idx);

  if (dim->arrsize < dim->numelems + 1) {
    if ((tmp = realloc (dim->arr, sizeof (double) * (dim->arrsize + 10))) == NULL)
      return (-1);
    dim->arr = tmp;
    dim->arrsize += 10;
  }
  if (idx < dim->numelems)
    memmove (dim->arr+idx+1, dim->arr+idx, sizeof (double) * (dim->numelems - idx));
  dim->arr[idx] = val;
  dim->numelems++;
  return (idx);
}


int find (st_dim *dim, double val, int *found)
{
  int va;

  *found = 0;
  if (dim->numelems == 0)
    return (dim->lastidx = 0);
  if (val == dim->arr[dim->lastidx]) {
    *found = 1;
    return (dim->lastidx);
  }

  if (val > dim->arr[dim->lastidx]) {
    for (va = dim->lastidx + 1; va < dim->numelems; va++) {
      if (val < dim->arr[va]) {
	return (dim->lastidx = va);
      } else if (val == dim->arr[va]) {
	*found = 1;
	return (dim->lastidx = va);
      }
    }
    return (dim->lastidx = dim->numelems);
    
  } else {   /*  val < dim->arr[dim->lastidx]  */
    for (va = 0; va <= dim->lastidx; va++) {
      if (val < dim->arr[va]) {
	return (dim->lastidx = va);
      } else if (val == dim->arr[va]) {
	*found = 1;
	return (dim->lastidx = va);
      }
    }
  }
  /* Should be impossible to reach this point */
  return (-1);
}


void multi_free (st_multidim *md)
{
  int va;

  for (va = 0; va < md->numdims; va++)
    free (md->dims[va].arr);
  free (md->dims);
  return;
}


void multi_dump (st_multidim *md)
{
  int va, vb;

  for (va = 0; va < md->numdims; va++) {
    printf ("dim %d:", va);
    for (vb = 0; vb < md->dims[va].arrsize; vb++) {
      if (vb < md->dims[va].numelems) {
	printf (" %5.2f", md->dims[va].arr[vb]);
      } else {
	printf (" .....");
      }
    }
    printf (" (lastidx %d, dimsize %d)\n", md->dims[va].lastidx, md->dims[va].dimsize);
  }
  printf ("totalelems %d\n", md->totalelems);
  return;
}
