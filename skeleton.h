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

  @Note The version information is for the whole file since it's 
  from svn.

  @Note There is no prototype for bFunction because it is only
  used within the module.

*/
#include <stdio.h>

int externallyDefinedVariable;   ///< Comment on extern (follows definition)

/**

  Brief description of public function ending with a period.

  Detailed multiline function description thru to the end
  of the comment. Emacs can fill indented paragraphs fine.

  @author Bill Valentine-Cooper - overall content.
  @author J. Random Coder - tons of corrections.
  @par Reviewer
     K. Random Reviewer on 2009-08-11.
     The Whole Group on 2009-08-12.

  @par Global Inputs

  Here we would discuss the globals that are referenced without
  modification by this function. We don't mention explicit parameters
  because they're picked-out by doxygen itself from the function
  definition.

  Note that multiple paragraphs in a section like this are fine.

  @par Global Outputs

  Here we would discuss any globals that are modified by this function.

  @return the new value of secondArg, unless firstArg wasn't positive.
  @retval -1 firstArg was less than zero
  @retval 0 firstArg was zero

  @sa relatedCode.c relatedHeader.h

*/
void aFunction (int firstArg,    ///< Some integer to be copied to the second argument.
               int &secondArg,   ///< The address of an integer to receive the first argument.
	       );

#endif
