#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <float.h>
#include <errno.h>

#define DPRIME_COL 1
#define THETA_COL  2
#define LR_COL     3

typedef struct {
  char name1[32],
    name2[32];
  int num,
    numcols,
    *datacols;
} st_marker;

typedef struct {
  double dprime,
    theta,
    lr;
} st_data;

typedef struct {
  int arrsize,
    numelems,
    lastidx;
  double *arr;
} st_ctlblock;

/* Default values that can be overridden on the command line */
int sexspecific = 0;   /* -s sex-specific thetas */
double prior = 0.02;   /* -p prior probability */
double weight = 0.95;  /* -w weighting factor */
double cutoff = 0.05;  /* -c theta cutoff */

/* Globally available for error messages */
int lineno = 0;

int get_marker_line (st_marker *marker, FILE *fp);
int get_header_line (st_marker *marker, FILE *fp);
int get_data_line (st_marker *marker, st_data *data, FILE *fp);
double calc_ppl (st_ctlblock *dprime, st_ctlblock *theta, double **lr);
double calc_ldppl (st_ctlblock *dprime, st_ctlblock *theta, double **lr);
double calc_ppld (st_ctlblock *dprime, st_ctlblock *theta, double **lr);
int parse_command_line (int argc, char **argv);
int insert (st_ctlblock *blk, double val);
int find (st_ctlblock *blk, double val);


main (int argc, char **argv)
{
  int argidx, ret, datalines, va, vb;
  char buff[256];
  FILE *fp;
  st_marker marker;
  st_data data;
  st_ctlblock dprime_ctl, theta_ctl;
  double **lr;

  argidx = parse_command_line (argc, argv);
  if (argidx >= argc) {
    fprintf (stderr, "missing file name\n");
    exit (-1);
  }
  if ((fp = fopen (argv[argidx], "r")) == NULL) {
    fprintf (stderr, "open '%s' failed, %s\n", argv[argidx], strerror (errno));
    exit (-1);
  }
  
  memset (&marker, 0, sizeof (st_marker));
  memset (&dprime_ctl, 0, sizeof (st_ctlblock));
  memset (&theta_ctl, 0, sizeof (st_ctlblock));

  if ((ret = get_marker_line (&marker, fp)) == 0) {
    fprintf (stderr, "file '%s' is empty\n", argv[argidx]);
    exit (-1);
  }
  if ((ret = get_header_line (&marker, fp)) == 0) {
    fprintf (stderr, "file '%s' ends unexpectedly at line %d\n", argv[argidx], lineno);
    exit (-1);
  }
  while ((ret = get_data_line (&marker, &data, fp)) == 1) {
    if (ret == 0) {
      fprintf (stderr, "file '%s' ends unexpectedly at line %d\n", argv[argidx], lineno);
      exit (-1);
    }
    if (insert (&dprime_ctl, data.dprime) == -1) {
      fprintf (stderr, "insert dprime %lf failed, %s\n", data.dprime, strerror (errno));
      exit (-1);
    }
    if (insert (&theta_ctl, data.theta) == -1) {
      fprintf (stderr, "insert theta failed, %s\n", strerror (errno));
      exit (-1);
    }
  }
  if (fseek (fp, 0, SEEK_SET) == -1) {
    fprintf (stderr, "fseek failed, %s\n", strerror (errno));
    exit (-1);
  }

  if ((lr = malloc (sizeof (double *) * dprime_ctl.numelems)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  for (va = 0; va < dprime_ctl.numelems; va++) {
    if ((lr[va] = malloc (sizeof (double) * theta_ctl.numelems)) == NULL) {
      fprintf (stderr, "malloc failed, %s\n", strerror (errno));
      exit (-1);
    }
  }

  printf ("%34s %6s %6s %6s\n", " ", "PPL","LD-PPL", "PPLD");

  while (ret = get_marker_line (&marker, fp) != 0) {
    if (get_header_line (&marker, fp) == 0) {
      fprintf (stderr, "file '%s' ends unexpectedly at line %d\n", argv[argidx], lineno);
      exit (-1);
    }
    datalines = 0;
    for (va = 0; va < dprime_ctl.numelems; va++) 
      for (vb = 0; vb < theta_ctl.numelems; vb++)
	lr[va][vb] = DBL_MAX;
    while (get_data_line (&marker, &data, fp) == 1) {
      if (((va = find (&dprime_ctl, data.dprime)) == -1) || (va >= dprime_ctl.numelems)) {
	fprintf (stderr, "unexpected dprime %lf at line %d\n", data.dprime, lineno);
	exit (-1);
      }
      if (((vb = find (&theta_ctl, data.theta)) == -1) || (vb >= theta_ctl.numelems)) {
	fprintf (stderr, "unexpected theta %lf at line %d\n", data.theta, lineno);
	exit (-1);
      }
      if (lr[va][vb] != DBL_MAX) {
	fprintf (stderr, "duplicate dprime %lf/theta %lf at line %d\n", data.dprime,
		 data.theta, lineno);
	exit (-1);
      }
      lr[va][vb] = data.lr;
      datalines++;
    }
    if (datalines != (va = dprime_ctl.numelems * theta_ctl.numelems)) {
      fprintf (stderr, "expected %d data lines, found %d, ending at line %d\n", va, 
	       datalines, lineno);
      exit (-1);
    }
    printf ("%2d %15s %15s %6.4f %.4f %.4f\n", marker.num, marker.name1, marker.name2, 
	    calc_ppl (&dprime_ctl, &theta_ctl, lr), calc_ldppl (&dprime_ctl, &theta_ctl, lr),
	    calc_ppld (&dprime_ctl, &theta_ctl, lr));
  }

  for (va = 0; va < dprime_ctl.numelems; va++) {
    free (lr[va]);
  }
  free (lr);
  free (dprime_ctl.arr);
  free (theta_ctl.arr);
  free (marker.datacols);
  exit (0);
}


double calc_ppl (st_ctlblock *dprime, st_ctlblock *theta, double **lr)
{
  int va, vb;
  double pre_int=0, post_int=0, integral, ppl, cutlr;

  if (((va = find (dprime, 0.0)) == -1) || (va >= dprime->numelems)) {
    fprintf (stderr, "can't calculate PPL, no dprime == 0\n");
    exit (-1);
  }

  for (vb = 1; vb < theta->numelems; vb++) {
    if (theta->arr[vb] <= cutoff) {
      pre_int += (lr[va][vb-1] + lr[va][vb]) * (theta->arr[vb] - theta->arr[vb-1]);
      
    } else if (theta->arr[vb-1] >= cutoff) {
      post_int += (lr[va][vb-1] + lr[va][vb]) * (theta->arr[vb] - theta->arr[vb-1]);
      
    } else {
      cutlr = lr[va][vb-1] +
	((cutoff - theta->arr[vb-1]) / (theta->arr[vb] - theta->arr[vb-1])) *
	(lr[va][vb] - lr[va][vb-1]);
      pre_int += (lr[va][vb-1] + cutlr) * (cutoff - theta->arr[vb-1]);
      post_int += (cutlr + lr[va][vb]) * (theta->arr[vb] - cutoff);
    }
  }

  /* The 0.5 here is factored out of the calculation of area. */
  pre_int *= (weight / cutoff) * 0.5;
  post_int *= ((1 - weight) / (0.5 - cutoff)) * 0.5;
  integral = pre_int + post_int;
  ppl = (prior * integral) / (prior * integral + (1 - prior));
  return (ppl);
}


double calc_ldppl (st_ctlblock *dprime, st_ctlblock *theta, double **lr)
{
  int va, vb;
  double pre_int=0, post_int=0, integral, ldppl, cutlr, l, r=0;

  for (va = 0; va < dprime->numelems; va++)
    r += lr[va][0];
  r /= dprime->numelems;

  for (vb = 1; vb < theta->numelems; vb++) {
    l = r;
    r = 0;
    for (va = 0; va < dprime->numelems; va++)
      r += lr[va][vb];
    r /= dprime->numelems;
    
    if (theta->arr[vb] <= cutoff) {
      pre_int += (l + r) * (theta->arr[vb] - theta->arr[vb-1]);
      
    } else if (theta->arr[vb-1] >= cutoff) {
      post_int += (l + r) * (theta->arr[vb] - theta->arr[vb-1]);
      
    } else {
      cutlr = l +
	((cutoff - theta->arr[vb-1]) / (theta->arr[vb] - theta->arr[vb-1])) * (r - l);
      pre_int += (l + cutlr) * (cutoff - theta->arr[vb-1]);
      post_int += (cutlr + r) * (theta->arr[vb] - cutoff);
    }
  }

  /* The 0.5 here is factored out of the calculation of area. */
  pre_int *= (weight / cutoff) * 0.5;
  post_int *= ((1 - weight) / (0.5 - cutoff)) * 0.5;
  integral = pre_int + post_int;
  ldppl = (prior * integral) / (prior * integral + (1 - prior));
  return (ldppl);
}


double calc_ppld (st_ctlblock *dprime, st_ctlblock *theta, double **lr)
{
  int va, vb, ndprimes = 0;
  double ld_pre=0, ld_post=0, le_pre=0, le_post=0, ppld, cutlr;
 
  for (va = 0; va < dprime->numelems; va++) {
    if (dprime->arr[va] == 0.0) {
      for (vb = 1; vb < theta->numelems; vb++) {
	if (theta->arr[vb] <= cutoff) {
	  le_pre += (lr[va][vb-1] + lr[va][vb]) * (theta->arr[vb] - theta->arr[vb-1]) *
	    (1 - prior);
	} else if (theta->arr[vb-1] >= cutoff) {
	  le_post += (lr[va][vb-1] + lr[va][vb]) * (theta->arr[vb] - theta->arr[vb-1]);
	} else {
	  cutlr = lr[va][vb-1] +
	    ((cutoff - theta->arr[vb-1]) / (theta->arr[vb] - theta->arr[vb-1])) *
	    (lr[va][vb] - lr[va][vb-1]);
	  le_pre += (lr[va][vb-1] + cutlr) * (cutoff - theta->arr[vb-1]) * (1 - prior);
	  le_post += (cutlr + lr[va][vb]) * (theta->arr[vb] - cutoff);
	}
      }
    } else {
      ndprimes++;
      for (vb = 1; vb < theta->numelems; vb++) {
	if (theta->arr[vb] <= cutoff) {
	  ld_pre += (lr[va][vb-1] + lr[va][vb]) * (theta->arr[vb] - theta->arr[vb-1]) * prior;
	} else if (theta->arr[vb-1] >= cutoff) {
	  /*ld_post += (lr[va][vb-1] + lr[va][vb]) * (theta->arr[vb] - theta->arr[vb-1]);*/
	} else {
	  cutlr = lr[va][vb-1] +
	    ((cutoff - theta->arr[vb-1]) / (theta->arr[vb] - theta->arr[vb-1])) *
	    (lr[va][vb] - lr[va][vb-1]);
	  ld_pre += (lr[va][vb-1] + cutlr) * (cutoff - theta->arr[vb-1]) * prior;
	  /*ld_post += (cutlr + lr[va][vb]) * (theta->arr[vb] - cutoff);*/
	}
      }
    }
  }

  /* The 0.5 here is factored out of the calculation of area */
  ld_pre *= (weight / cutoff) * 0.5;
  ld_pre /= ndprimes;
  /*ld_post *= ((1 - weight) / (0.5 - cutoff)) * 0.5;*/
  /*ld_post /= ndprimes;*/
  le_pre *= (weight / cutoff) * 0.5;
  le_post *= ((1 - weight) / (0.5 - cutoff)) * 0.5;
  ppld = (ld_pre + ld_post) / (ld_pre + ld_post + le_pre + le_post);
  return (ppld);
}


int get_marker_line (st_marker *marker, FILE *fp)
{
  char buff[256];
  char *pa, *pb;

  if (fgets (buff, 256, fp) == NULL) {
    if (feof (fp))
      return (0);
    fprintf (stderr, "read error at line %d, %s\n", lineno, strerror (errno));
    exit (-1);
  }
  lineno++;

  if (((pa = strtok_r (buff, " \t\n", &pb)) == NULL) ||
      (strcmp (pa, "#") != 0))
    return (-1);
  
  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL)
    return (-1);
  if (((marker->num = (int) strtol (pa, NULL, 10)) == 0) && (errno != 0))
    return (-1);

  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL)
    return (-1);
  strcpy (marker->name1, pa);

  if ((pa = strtok_r (NULL, " \t\n", &pb)) == NULL)
    return (-1);
  strcpy (marker->name2, pa);

  return (1);
}


int get_header_line (st_marker *marker, FILE *fp)
{ 
  char buff[256];

  if (fgets (buff, 256, fp) == NULL) {
    if (feof (fp))
      return (0);
    fprintf (stderr, "read error at line %d, %s\n", lineno, strerror (errno));
    exit (-1);
  }
  lineno++;

  /* We'll only parse one header line per file, and require that the
   * format be consistant within the file. If this pointer is non-null,
   * we've parsed a header previsously, so just return.
   */
  if (marker->datacols != NULL)
    return (1);

  /* default values for expected input format. worry about dynamically 
   * determining input format some other day.
   */
  if ((marker->datacols = malloc (sizeof (int) * 5)) == NULL) {
    fprintf (stderr, "malloc failed, %s\n", strerror (errno));
    exit (-1);
  }
  memset (marker->datacols, 0, sizeof (int) * 5);
  marker->numcols = 5;
  marker->datacols[0] = DPRIME_COL;
  marker->datacols[1] = THETA_COL;
  marker->datacols[4] = LR_COL;

  return (1);
}


int get_data_line (st_marker *marker, st_data *data, FILE *fp)
{
  char buff[256];
  char *pa, *pb;
  long start;
  double da;
  int va = 0;

  if ((start = ftell (fp)) == -1) {
    fprintf (stderr, "ftell failed, %s\n", strerror (errno));
    exit (-1);
  }
  if (fgets (buff, 256, fp) == NULL) {
    if (feof (fp))
      return (0);
    fprintf (stderr, "read error at line %d, %s\n", lineno, strerror (errno));
    exit (-1);
  }
  lineno++;
  
  if ((pa = strtok_r (buff, " (),\t\n", &pb)) == NULL) {
    fprintf (stderr, "unexpectedly short line %d\n", lineno);
    exit (-1);
  }

  /* If it's a new marker line, back the file pointer up, and return a 
   * non-error, non-dataline value.
   */
  if (strcmp (pa, "#") == 0) {
    if (fseek (fp, start, SEEK_SET) == -1) {
      fprintf (stderr, "fseek failed, %s\n", strerror (errno));
      exit (-1);
    }
    return (2);
  }

  while (1) {
    if (marker->datacols[va] == DPRIME_COL) {
      data->dprime = strtod (pa, NULL);
    } else if (marker->datacols[va] == THETA_COL) {
      data->theta = strtod (pa, NULL);
    } else if (marker->datacols[va] == LR_COL) {
      data->lr = strtod (pa, NULL);
    }
    if (++va >= marker->numcols)
      break;
    if ((pa = strtok_r (NULL, " (),\t\n", &pb)) == NULL) {
      fprintf (stderr, "unexpectedly short line %d\n", lineno);
      exit (-1);
    }
  }

  return (1);
}


int parse_command_line (int argc, char **argv)
{
  int arg;

  opterr = 0;
  while ((arg = getopt (argc, argv, "sp:w:c:")) != -1) {
    if (arg == 's') {
      sexspecific = 1;
    } else if (arg == 'p') {
      prior = strtod (optarg, NULL);
    } else if (arg == 'w') {
      weight = strtod (optarg, NULL);
    } else if (arg == 'c') {
      cutoff = strtod (optarg, NULL);
    } else {
      fprintf (stderr, "bad option\n");
      exit (-1);
    }
  }
  return (optind);
}


int insert (st_ctlblock *blk, double val)
{
  double *tmp;
  int idx;
  
  if (blk->arrsize == 0) {
    if ((blk->arr = malloc (sizeof (double) * 10)) == NULL) 
      return (-1);
    blk->arrsize = 10;
    blk->numelems = 1;
    blk->lastidx = 0;
    blk->arr[0] = val;
    return (0);
  } else {
    if ((idx = find (blk, val)) == -1)
      return (-1);
    if (idx < blk->numelems)
      return (idx);
    if (idx >= blk->arrsize) {
      if ((tmp = realloc (blk->arr, sizeof (double) * (blk->arrsize + 10))) == NULL)
	return (-1);
      blk->arr = tmp;
      blk->arrsize += 10;
    }
    blk->arr[idx] = val;
    blk->numelems++;
    return (blk->lastidx = idx);
  }
}


int find (st_ctlblock *blk, double val)
{
  int va;
  
  if (blk->numelems == 0)
    return (-1);
  if (val == blk->arr[blk->lastidx])
    return (blk->lastidx);
  
  if (val > blk->arr[blk->lastidx]) {
    for (va = blk->lastidx + 1; va < blk->numelems; va++) {
      if (val == blk->arr[va])
	return (blk->lastidx = va);
    }
    if (val > blk->arr[blk->numelems - 1])
      return (blk->numelems);
    return (-1);
    
  } else {   /*  val < blk->arr[blk->lastidx]  */
    for (va = 0; va < blk->lastidx; va++) {
      if (val == blk->arr[va])
	return (blk->lastidx = va);
    }
    return (-1);
  }
}
