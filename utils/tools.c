
/**********************************************************************
 * Copyright 2008, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

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
