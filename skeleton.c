/**
@file skeleton.c

  Brief internal-use-only module description ending with a period.

  Detailed internal-use-only module description.
  This will be treated like HTML, so embed any format tags you might
  need.  The exception is paragraph breaks, which are not needed.

  Copyright 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

  Note that the version information is for the whole file since it's 
  from svn.

*/
#include "skeleton.h"

int moduleGlobalVariable,       ///< Comment on global
  anotherGlobal,        ///< Another comment
  lastOne;      ///< Comment on last global

/*

  Documented in header.

*/
void aFunction (int firstArg, int &secondArg)
{
  /*
    This is not really intended to be an example
    of how to write code. This comment will not
    show up in the documentation.
  */
  if (firstArg < 0)
    return -1;
  else  // Neither will this comment.
    if (firstArg == 0)
      return 0;

  /**
     @par Code Comment.

     This is a section discussing the algorithm used
     way down in the code somewhere. You don't have to
     use "Code Comment" up there as the paragraph name.
     I just chose that to make it stand-out from the
     rest of the documentation which came from the 
     header. Otherwise, it might seem confusingly 
     out-of-context.

  */
  *secondArg = firstArg;
  return secondArg;
}


/**

  Brief function description ending with a period.

  Detailed multiline function description thru to the end
  of the comment. Emacs can fill indented paragraphs fine.

  @author Bill Valentine-Cooper - overall content.
  @author J. Random Coder - tons of corrections.
  @par Reviewer
     K. Random Reviewer on 2009-08-11.

  @par Global Inputs

  Here we would discuss the globals that are referenced without
  modification by this function. We don't mention explicit parameters
  because they're picked-out by doxygen itself from the function
  definition.

  Note that multiple paragraphs in a section like this are fine.

  @par Global Outputs

  Here we would discuss any globals that are modified by this function.

  @return nothing

*/
void
bFunction () {
  return;
}
