#ifndef __SKELETON_H__
#define __SKELETON_H__

/**
@file skeleton.h

  @see skeleton.c
  
  Copyright &copy; 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

  @note There is no prototype for bFunction because it is private,
  i.e. only used within the module.

*/
#include <stdio.h>

/// Comment on extern (preceeds target)
int externallyDefinedVariable;

void aFunction (int firstArg,    ///< This comment will be overridden by equivalent in function source.
               int &secondArg   ///< The address of an integer to receive the first argument in code.
	       );

#endif
