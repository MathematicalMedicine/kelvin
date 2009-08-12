#ifndef SKELETON_H
#define SKELETON_H

/**
@file skeleton.h

  Brief public module description ending with a period.

  Detailed public module description. This will be treated like HTML, 
  so embed any format tags you might
  need.  The exception is paragraph breaks, which are not needed.

  Note that the detailed overall module description as well as the 
  documentation for global functions and variables is in the module header,
  because if we were writing a library that we distributed in binary
  form, we would provide no source, but only headers.

  Copyright 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

  @note The version information is for the whole file since it's 
  from svn.

  @note There is no prototype for bFunction because it is only
  used within the module.

*/
#include <stdio.h>

/// Comment on extern (preceeds target)
int externallyDefinedVariable;

void aFunction (int firstArg,    ///< This comment will be overridden by equivalent in function source.
               int &secondArg   ///< The address of an integer to receive the first argument in code.
	       );

#endif
