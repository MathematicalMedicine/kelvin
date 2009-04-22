
/**********************************************************************
 * Kelvin utilities.
 * Alberto Maria Segre
 *
 * Copyright 2005, The University of Iowa.  All rights reserved.
 * Permission is hereby given to use and reproduce this software 
 * for non-profit educational purposes only.
 **********************************************************************/
#include "utils.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

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

  logFlag[LOGERROR] = -1; // Give errors for everything by default
  logFlag[LOGWARNING] = -1; // Give warnings for everything by default
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
