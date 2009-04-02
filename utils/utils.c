/**********************************************************************
 * Kelvin utilities.
 * Alberto Maria Segre, et al.
 *
 * Copyright 2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "utils.h"

/**********************************************************************
 * Kelvin log facility. Keep an array of bitvectors, where each
 * bitvector corresponds to a level of output verbosity and each bit
 * corresponds to an independent "type" or stream of output.
 *
 * Levels are defined from 0 (most serious, inevitably fatal, and
 * always produced) to MAXLOGLEVELS (least serious). Commonly used
 * names for the levels are defined in log.h.
 *
 * Since level 0 errors are always produced, level k error bits are
 * actually set in logFlag[k-1] (in other words, we don't need flag
 * bits for level 0, so level 1 flag bits reside in logFlag[0]
 * instead).
 *
 * Types of output are also defined in log.h and are specific to
 * Kelvin. Note that with 32 bit ints you are limited to 32 different
 * output types; if you need more, you'll have to redefine logFlag as
 * something longer than 32 bits.
 **********************************************************************/
unsigned int logFlag[MAXLOGLEVELS];

/* Initialize the log facility. LOGDEFAULT errors are on by
 * default. */
void
logInit ()
{
  int i;

  for (i = 0; i < MAXLOGLEVELS; i++)
    logFlag[i] = LOGDEFAULT;
}

/* Set a particular type of output to a given level severity or
 * greater (which, internally, means this output will be produced if
 * level index is actually less than the specified level minus 1,
 * since (a) 0 is the most severe error, and we use "off by 1"
 * indexing to account for the fact that level 0 errors are always
 * produced). */
void
logSet (unsigned int type, int level)
{
  int i;

  for (i = 0; i < level; i++)
    logFlag[i] |= type;
}

/* Print an output message if the specified type of output is turned
 * on at this level of output. Careful of off by 1 indexing of levels
 * here! */
void
logMsg (unsigned int type, int level, const char *format, ...)
{
  va_list argp;

  /* Initialize the variable arguments list. */
  va_start (argp, format);
  /* Check to see if the criteria for making the message appear are
   * met. Level 0 errors are always produced. */
  if ((level == 0) || (type & logFlag[level - 1]))
    vfprintf (stderr, format, argp);
  /* Close the variable arguments list. */
  va_end (argp);

  /* If this was a fatal error, dump. Level 0 errors are always checked
   * and are always fatal. */
  if (level == LOGFATAL)
    exit (ERROR);
}


/* Wrappers for malloc() and friends, with built in error checking
 * and logging.
 */
void *
MALLOC (char *description, size_t size)
{
  void *ptr = NULL;
  ptr = malloc (size);
  if (ptr == NULL)
    {
      fprintf (stderr, "Memory allocation failed.\n");
    }
  KASSERT (ptr != NULL, "Can't allocate memory for %s of size %d.\n",
	   description, size);
  KLOG (LOGMEMORY, LOGDEBUG, "Allocate %p of size %d - %s \n",
	ptr, size, description);
  return ptr;
}

void
FREE (char *description, void *ptr)
{
  KLOG (LOGMEMORY, LOGDEBUG, "Free %p - %s \n", ptr, description);
  free (ptr);
}

void *
REALLOC (char *description, void *ptr, size_t size)
{
  void *oldPtr;
  void *newPtr = NULL;

  oldPtr = ptr;
  newPtr = realloc (ptr, size);
  if (newPtr == NULL)
    {
      fprintf (stderr, "Memory reallocation failed.\n");
    }
  if (ptr == NULL)
    {
      memset (newPtr, 0, size);
    }

  KASSERT (newPtr != NULL, "Can't re-allocate memory for %s of size %d.\n",
	   description, size);
  KLOG (LOGMEMORY, LOGDEBUG, "Reallocate %p(from %p) of size %d - %s \n",
	newPtr, oldPtr, size, description);
  return newPtr;
}


/* Input handling/checking routines, used when reading configuration
 * and data files.
 */
int
is_line_blank (char *line)
{
  int pos = 0;			/* current position in the line */

  while (line[pos] != '\0' &&
	 (line[pos] == ' ' || line[pos] == '\t' || line[pos] == '\n'))
    pos++;

  /* if we have reached the end of the line, then it's a blank line */
  if (pos == strlen (line))
    return 1;
  else
    return 0;
}

int
is_line_comment (char *line)
{
  int pos = 0;			/* current position in the line */

  while (line[pos] != '\0' &&
	 (line[pos] == ' ' || line[pos] == '\t' || line[pos] == '\n'))
    pos++;

  /* if we have reached the end of the line, then it's a blank line */
  if (pos == strlen (line))
    return 0;
  else if (line[pos] == '#')
    /* we got a comment line */
    return 1;
  else
    return 0;
}

int
is_line_blank_or_comment (char *line)
{
  int pos = 0;			/* current position in the line */

  while (line[pos] != '\0' &&
	 (line[pos] == ' ' || line[pos] == '\t' || line[pos] == '\n'))
    pos++;

  /* if we have reached the end of the line, then it's a blank line */
  if (pos == strlen (line))
    return 1;
  if (line[pos] == '#')
    /* we got a comment line */
    return 1;

  /* not a comment or blank line */
  return 0;
}

char *
get_nonblank_line (char *pLine, int maxLen, FILE * fp, int *pLineNo)
{
  while (fgets (pLine, maxLen, fp))
    {
      if (pLineNo != NULL)
	(*pLineNo)++;
      if (is_line_blank_or_comment (pLine))
	continue;
      return pLine;
    }
  /* reaching end of file */
  return NULL;
}

#define BUFF_INCR 512

char *
fgetlongs (char **buff, int *bufflen, FILE * fp)
{
  int offset = 0;

  do
    {
      if (*bufflen == 0)
	{
	  if ((*buff = (char *) malloc (BUFF_INCR)) == NULL)
	    return (NULL);
	  *bufflen = BUFF_INCR;
	}
      else if (offset != 0)
	{
	  if ((*buff =
	       (char *) realloc (*buff, *bufflen + BUFF_INCR)) == NULL)
	    return (NULL);
	  *bufflen += BUFF_INCR;
	}
      if (fgets (*buff + offset, *bufflen - offset, fp) == NULL)
	{
	  if (offset != 0)
	    return (*buff);
	  return (NULL);
	}
    }
  while (((offset = strlen (*buff)) == (*bufflen - 1))
	 && ((*buff)[offset - 1] != 0x0a));

  return (*buff);
}

/* Symbols that are only useful to permuteLine(), so no need to polute the
 * namespace by putting them in the header file.
 */
#define STARTOFLINE  0
#define INSTRING     1
#define INWHITESPACE 2
#define INSEPARATOR  3
#define INCOMPARATOR 4

/* Permutes a line of input. Leading and trailing whitespace is deleted;
 * whitespace that brackets 'separator' characters (',', '-', ';', ':') is
 * deleted; whitespace between 'comparator' characters ('!', '>'. '=') is deleted,
 * all other whitespace is reduced to a single space; groups of comparator
 * characters are forcibly bracketed by single spaces; comment characters
 * (pound sign) and all subsequent characters up to end of line are deleted.
 * 
 */
int permuteLine (char *line, int maxlength)
{
  int state, compcount=0, va, vb;

  va = vb = 0;
  state = STARTOFLINE;

  /* First loop: strip whitespace around separators AND comparators; reduce whitespace 
   * around strings to a single space. At the same time, count up how many spaces we'll
   * need to put around groups of comparators in the second loop.
   */
  while (1) {
    //printf ("va %d is '%c', vb %d is '%c' state is %d -> ", va, line[va], vb, line[vb], state);
    if (index (" \t", line[vb]) != NULL) {
      if (state == INSTRING) 
	line[va++] = ' ';
      else if (state == INCOMPARATOR)
	compcount++;
      vb++;
      if (state != STARTOFLINE)
	state = INWHITESPACE;
      
    } else if (index ("-:;,", line[vb]) != NULL) {
      if (state == INWHITESPACE)
	line[va-1] = line[vb++];
      else 
	line[va++] = line[vb++];
      if (state == INCOMPARATOR)
	compcount++;
      state = INSEPARATOR;

    } else if (line[vb] == '#') {
      if (state == INWHITESPACE)
	va--;
      line[va] = '\0';
      break;

    } else if ((index ("\n\r", line[vb]) != NULL) || (line[vb] == '\0')) {
      if (state == INWHITESPACE)
	va--;
      line[va] = '\0';
      break;

    } else if (index ("=>!", line[vb]) != NULL) {
      if (state == INWHITESPACE)
	line[va-1] = line[vb++];
      else 
	line[va++] = line[vb++];
      if (state != INCOMPARATOR)
	compcount++;
      state = INCOMPARATOR;
      
    } else {
      line[va++] = line[vb++];
      if (state == INCOMPARATOR)
	compcount++;
      state = INSTRING;
    }
    //printf ("va %d is '%c', vb %d is '%c' state is %d, \n", va, line[va], vb, line[vb], state);
  }
  //printf ("all done: va %d, vb %d, compcount %d, line '%s'\n", va, vb, compcount, line);

  if (compcount == 0)
    return (0);
  if (strlen (line) + compcount + 1 > maxlength)
    return (-1);
    
  /* Second loop: insert spaces around groups of comparators. We simplify the states;
   * now strings are "everything that's not a comparator".
   */
  va = strlen (line);
  vb = va + compcount;
  line[va--] = ' ';
  line[vb--] = '\0';
  state = INSTRING;

  while (va != vb) {
  printf ("va %2d, vb %2d, state %d, line '%s' -> ", va, vb, state, line);
    if (index ("=>!", line[va]) != NULL) {
      if (state == INSTRING)
	line[vb--] = ' ';
      else 
	line[vb--] = line[va--];
      state = INCOMPARATOR;

    } else {
      if (state == INCOMPARATOR)
	line[vb--] = ' ';
      else
	line[vb--] = line[va--];
      state = INSTRING;
    }
    printf ("'%s'\n", line);
  }  

  return (0);
}
