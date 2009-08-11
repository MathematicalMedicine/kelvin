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

@author Bill Valentine-Cooper.

@param[in] firstArg - some integer to be copied to the second argument.
@param[out] secondArg - the address of an integer to receive the first argument.
@return nothing.

*/
void function (int firstArg,    ///< Description of argument
               int &secondArg,   ///< Description of next argument
	       )
{
  *secondArg = firstArg;
  return;
}
