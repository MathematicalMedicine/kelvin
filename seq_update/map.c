#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <map.h>

st_mapchr *mapchrs = NULL;
int nummapchrs = 0;
int nextchridx=0, nextmarkeridx=0;
char *mapfile;

void insert_mapmarker (st_mapmarker *mrk);
int mapchr_compare (const void *p1, const void *p2);
char *strlower (char *str);


void read_map (char *filename)
{
  int numcols=0, *datacols=NULL, lineno=0, colno;
  char buff[BUFFLEN], *token, *pb, *endptr;
  FILE *fp;
  st_mapmarker mrk;
  
  mapfile = filename;
  if ((fp = fopen (mapfile, "r")) == NULL) {
    fprintf (stderr, "open '%s' failed, %s\n", mapfile, strerror (errno));
    exit (-1);
  }

  /* I don't care what your stinky map function is */
  while (1) {
    if (fgets (buff, BUFFLEN, fp) == NULL) {
      if (feof (fp))
	fprintf (stderr, "mapfile '%s' ends unexpectedly at line %d\n", mapfile, lineno);
      else 
	fprintf (stderr, "error reading '%s' at line %d, %s\n", mapfile, lineno,
		 strerror (errno));
      exit (-1);
    }
    lineno++;
    if (strstr (strlower (buff), "mapfunction") == NULL)
      break;
  }
  
  /* We'll try to be flexible about file format, so long as there are column headers that
   * allow us to identify the columns. 
   */

  token = strtok_r (buff, " \t\n", &pb);
  while (token != NULL) {
    numcols++;
    if ((datacols = (int *) realloc (datacols, sizeof (int) * numcols)) == NULL) {
      fprintf (stderr, "realloc datacols failed, %s\n", strerror (errno));
      exit (-1);
    }

    if (strncmp (token, "chr", 3) == 0) {
      /* The chromosome column */
      datacols[numcols - 1] = MAP_CHR_COL;
      
    } else if ((strncmp (token, "name", 4) == 0) || (strncmp (token, "marker", 6) == 0)) {
      /* The marker column */
      datacols[numcols - 1] = MAP_NAME_COL;
      
    } else if (strstr (token, "female") != NULL) {
      /* sex-specific position for ladies */
      datacols[numcols - 1] = MAP_FEMALEPOS_COL;

    } else if (strstr (token, "male") != NULL) {
      /* sex-specific position for gents */
      datacols[numcols - 1] = MAP_MALEPOS_COL;

    } else if ((strncmp (token, "basepair", 8) == 0) || (strstr (token, "phys") != NULL)) {
      /* physical or basepair position */
      datacols[numcols - 1] = MAP_BASEPAIR_COL;
      
    } else if ((strstr (token, "sex") != NULL) || (strstr (token, "ave") != NULL) ||
	       (strstr (token, "agv") != NULL) || (strstr (token, "pos") != NULL) ||
	       (strncmp (token, "kosambi", 7) == 0)) {
      /* The sex-averaged position column */
      datacols[numcols - 1] = MAP_AVGPOS_COL;
      
    } else {
      /* Something else */
      datacols[numcols - 1] = 0;
    }
    token = strtok_r (NULL, " \t\n", &pb);
  }

  while (fgets (buff, BUFFLEN, fp) != NULL) {
    lineno++;
    memset (&mrk, 0, sizeof (st_mapmarker));
    mrk.basepair = (long) (mrk.malepos = mrk.femalepos = -1.0);
    
    token = strtok_r (buff, " \t\n", &pb);
    for (colno = 0; colno < numcols; colno++) {
      if (token == NULL) {
	fprintf (stderr, "short line in '%s', line %d\n", mapfile, lineno);
	exit (-1);
      }
      switch (datacols[colno]) {
      case MAP_CHR_COL:
	if (strlen (token) > MAX_MAP_CHR_LEN-1) {
	  fprintf (stderr, "chromosome identifier '%s' in mapfile '%s' is too long\n", token,
		   mapfile);
	  exit (-1);
	}
	endptr = strcpy (mrk.chr, token);
	break;
      case MAP_NAME_COL:
	if (strlen (token) > MAX_MAP_NAME_LEN-1) {
	  fprintf (stderr, "marker name '%s' in mapfile '%s' is too long\n", token, mapfile);
	  exit (-1);
	}
	endptr = strcpy (mrk.name, token);
	break;
      case MAP_ALTNAME_COL:
	if (strlen (token) > MAX_MAP_NAME_LEN-1) {
	  fprintf (stderr, "marker name '%s' in mapfile '%s' is too long\n", token, mapfile);
	  exit (-1);
	}
	endptr = strcpy (mrk.altname, token);
	break;
      case MAP_AVGPOS_COL:
	mrk.avgpos = strtod (token, &endptr);
	break;
      case MAP_MALEPOS_COL:
	mrk.malepos = strtod (token, &endptr);
	break;
      case MAP_FEMALEPOS_COL:
	mrk.femalepos = strtod (token, &endptr);
	break;
      case MAP_BASEPAIR_COL:
	mrk.basepair = strtol (token, &endptr, 10);
	break;
      }
      if (token == endptr) {
	fprintf (stderr, "illegal data in line %d in '%s'\n", lineno, mapfile);
	exit (-1);
      }
      token = strtok_r (NULL, " \t\n", &pb);
    }
    if (token != NULL) {
      fprintf (stderr, "extra data in line %d in '%s'\n", lineno, mapfile);
      exit (-1);
    }
    insert_mapmarker (&mrk);
  }
  if (! feof (fp)) {
    fprintf (stderr, "error reading '%s' at line %d, %s\n", mapfile, lineno, strerror (errno));
    exit (-1);
  }
  if (nummapchrs > 1)
    qsort (mapchrs, nummapchrs, sizeof (st_mapchr), mapchr_compare);
  fclose (fp);
  free (datacols);
  return;
  
}


void insert_mapmarker (st_mapmarker *mrk)
{
  int va;
  st_mapchr *chr_p = NULL;

  for (va = 0; va < nummapchrs; va++) {
    if (strcmp (mapchrs[va].chr, mrk->chr) == 0) {
      chr_p = &mapchrs[va];
      break;
    }
  }
  if (chr_p == NULL) {
    if ((mapchrs = (st_mapchr *) realloc (mapchrs, sizeof (st_mapchr) * (nummapchrs + 1))) == NULL) {
      fprintf (stderr, "realloc mapchrs failed, %s\n", strerror (errno));
      exit (-1);
    }
    chr_p = &mapchrs[nummapchrs];
    nummapchrs++;
    memset (chr_p, 0, sizeof (st_mapchr));
    strcpy (chr_p->chr, mrk->chr);
  }

  if ((chr_p->nummarkers != 0) && (chr_p->markers[chr_p->nummarkers-1].avgpos > mrk->avgpos)) {
    fprintf (stderr, "position '%f' out of order in '%s'\n", mrk->avgpos, mapfile);
    exit (-1);
  }
  if ((chr_p->markers = realloc (chr_p->markers, sizeof (st_mapmarker) * (chr_p->nummarkers+1))) == NULL) {
    fprintf (stderr, "realloc markers failed, %s\n", strerror (errno));
    exit (-1);
  }
  memcpy (&chr_p->markers[chr_p->nummarkers], mrk, sizeof (st_mapmarker));
  chr_p->nummarkers++;
  return;
}


st_mapmarker *find_mapmarker (char *chr, char *name)
{
  int va, vb;
  st_mapchr *chr_p = NULL;

  for (va = 0; va < nummapchrs; va++) {
    if (strcmp (mapchrs[va].chr, chr) == 0) {
      chr_p = &mapchrs[vb = va];
      break;
    }
  }
  if (chr_p == NULL)
    return (NULL);
  for (va = 0; va < chr_p->nummarkers; va++) {
    if (strcmp (chr_p->markers[va].name, name) == 0) {
      nextchridx = vb;
      nextmarkeridx = va+1;
      return (&(chr_p->markers[va]));
    }
  }
  return (NULL);
}


st_mapmarker *next_mapmarker (st_mapmarker *mrk)
{
  if (nummapchrs == 0)
    return (NULL);

  if ((mrk != NULL) && (find_mapmarker (mrk->chr, mrk->name) == NULL))
    return (NULL);

  if ((nextchridx < nummapchrs) && (nextmarkeridx < mapchrs[nextchridx].nummarkers)) {
    return (&(mapchrs[nextchridx].markers[nextmarkeridx++]));
  } else if ((nextchridx + 1) < nummapchrs) {
    nextchridx++;
    nextmarkeridx = 0;
    return (&(mapchrs[nextchridx].markers[nextmarkeridx++]));
  } else {
    nextchridx = 0;
    nextmarkeridx = 0;
    return (NULL);
  }
}


void dump_map ()
{
  int va, vb;
  st_mapmarker *ptr;

  for (va = 0; va < nummapchrs; va++) {
    printf ("chr %s, %d markers\n", mapchrs[va].chr, mapchrs[va].nummarkers);
    for (vb = 0; vb < mapchrs[va].nummarkers; vb++) {
      ptr = &(mapchrs[va].markers[vb]);
      printf ("  %s pos %.6f (male %.6f/female %.6f) basepair %ld\n",
	      ptr->name, ptr->avgpos, ptr->malepos, ptr->femalepos, ptr->basepair);
    }
  }
  return;
}


void free_map ()
{
  int va;

  for (va = 0; va < nummapchrs; va++) {
    if (mapchrs[va].nummarkers > 0)
      free (mapchrs[va].markers);
  }
  if (nummapchrs > 0)
    free (mapchrs);
  return;
}


int mapchr_compare (const void *p1, const void *p2)
{
  return (strcmp (((st_mapchr *) p1)->chr, ((st_mapchr *) p2)->chr));
}


char *strlower (char *str)
{
  int va=0;

  while (str[va] != '\0') {
    str[va] = tolower (str[va]);
    va++;
  }
  return (str);
}
