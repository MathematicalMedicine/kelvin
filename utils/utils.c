/**
@file utils.c

  Utility stuff included by everything.

  String utilities, and otherwise useful macros.
  Most of it was originally written by Alberto Maria Segre.

  Copyright &copy; 2010, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> // For index on some platforms
#include <stdarg.h>
#include <limits.h>
#include <ctype.h>
#include "utils.h"

/* Input handling/checking routines, used when reading configuration
 * and data files.
 */
int
is_line_blank (char *line)
{
  unsigned int pos = 0;			/* current position in the line */

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
  unsigned int pos = 0;			/* current position in the line */

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
  unsigned int pos = 0;			/* current position in the line */

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

char *
strlower (char *str)
{
  int va=0;
  
  while (str[va] != '\0') {
    str[va] = tolower (str[va]);
    va++;
  }
  return (str);
}

#define BUFF_INCR 512UL

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

/* Symbols that are only useful to permuteLine(), so no need to pollute the
 * namespace by putting them in the header file.
 */
#define STARTOFLINE  0
#define INSTRING     1
#define INWHITESPACE 2
#define INSEPARATOR  3
#define INCOMPARATOR 4

/* Permutes a line of input:
 * - Leading and trailing whitespace is deleted
 * - Whitespace that brackets 'separator' characters ('-', ';', ':') is
 *   deleted, with the exeception of a '-' that is preceeded by a non-number
 *   and whitepace, and followed by a number
 * - Whitespace between 'comparator' characters (',', '!', '>'. '=') is deleted,
 *   and groups of comparator characters are forcibly bracketed by a single space
 * - All other whitespace is reduced to a single space
 * - Comment characters (pound sign) and all subsequent characters up to end
 *   of line are deleted.
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
    //printf ("%*sv (va %d)\n", va+1, " ", va);
    //printf (" %s", line);
    //printf ("%*s^ (vb %d, char %d, state %d, compcount %d)\n", vb+1, " ", vb, line[vb],
    //        state, compcount);

    /* It's important that we check first for a NULL character. index() will 
     * return success when searching any string for a NULL character, so we have
     * to check that case to have any confidence in index() later on.
     */
    if ((line[vb] == '\0') || (index ("\n\r", line[vb]) != NULL)) {
      if (state == INWHITESPACE)
	va--;
      line[va] = '\0';
      break;
      
    } else if (index (" \t", line[vb]) != NULL) {
      if (state == INSTRING) {
	line[va++] = ' ';
	if (state != STARTOFLINE)
	  state = INWHITESPACE;
      }
      vb++;
      
    } else if (line[vb] == '-') {  /* have to handle '-' separately */
      if ((state == INWHITESPACE) && (va > 1) && isdigit ((int) line[va-2])) {
	line[va-1] = line[vb++];
      } else 
	line[va++] = line[vb++];
      if (state == INCOMPARATOR)
	compcount++;
      state = INSEPARATOR;

    } else if (index (":;", line[vb]) != NULL) {   /* this used to include ',' */
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

    } else if (index ("=>!,", line[vb]) != NULL) {  /* ',' has been added here */
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
    //printf ("----\n");
  }
  //printf ("%*sv (va %d)\n", va+1, " ", va);
  //printf (" %s\n", line);
  //printf ("%*s^ (vb %d, state %d, compcount %d, final)\n", vb+1, " ", vb, state, compcount);
  //printf ("+++++\n");

  if (compcount == 0)
    return (0);
  if (((int) strlen (line)) + compcount + 1 > maxlength)
    return (-1);
    
  /* Second loop: insert spaces around groups of comparators. We simplify the states;
   * now "string" means "everything that's not a comparator".
   */
  va = strlen (line);
  vb = va + compcount;
  line[va--] = ' ';
  line[vb--] = '\0';
  state = INSTRING;

  while (va != vb) {
    //printf ("va %2d, vb %2d, state %d, line '%s' -> ", va, vb, state, line);
    if (index ("=>!,", line[va]) != NULL) {
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
    //printf ("'%s'\n", line);
  }  

  return (0);
}
