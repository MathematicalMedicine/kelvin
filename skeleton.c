/**
@file skeleton.c

  Brief module description ending with a period.

  Detailed multiline module description thru to the end of the comment.
  This will be treated like HTML, so embed any format tags you might
  need.  The exception is paragraph breaks, which are not needed.

  Be sure to include the actual module (or include) filename up at the
  top of the block.

  Copyright 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

  Note that the version information is for the whole file since it's 
  from svn.

*/
#include <standard header.h>
#include "local header.h"

extern int externallyDefinedVariable;   ///< Comment on extern (follows definition)

int moduleGlobalVariable,       ///< Comment on global
  anotherGlobal,        ///< Another comment
  lastOne;      ///< Comment on last global

/**

  Brief function description ending with a period.

  Detailed multiline function description thru to the end
  of the comment. Emacs can fill indented paragraphs fine.

  @author Bill Valentine-Cooper - overall content.
  @author J. Random Coder - tons of corrections.
  @par Reviewer K. Random Reviewer on 2009-08-12

  @par Inputs (global and explicit)

  Here we would discuss the globals that are referenced by this function as well
  as the explicit parameters.

  Multiple paragraphs are fine.

  @par Outputs (global and explicit)

  Here we would discuss any globals that are modified by this function, as well
  as any explicit parameters modified (output parameters).

  @return the new value of secondArg.
  @retval -1 firstArg was less than zero
  @retval 0 firstArg was zero

  @sa relatedCode.c relatedHeader.h

*/
void function (int firstArg,    ///< Some integer to be copied to the second argument.
               int &secondArg,   ///< The address of an integer to receive the first argument.
	       )
{
  if (firstArg < 0)
    return -1;
  else
    if (firstArg == 0)
      return 0;
  *secondArg = firstArg;
  return secondArg;
}
