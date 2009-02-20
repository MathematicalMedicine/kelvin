/**
@file skeleton.c

  Brief module description ending with a period.

  Detailed multiline module description thru to the end of the comment.
  This will be treated like HTML, so embed any format tags you might
  need.  The exception is paragraph breaks, which are not needed.

  Be sure to include the actual module (or include) filename up at the
  top of the block.

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

*/
void function (int firstArg,    ///< Description of argument
               int secondArg,   ///< Description of next argument
  )
{
  return;
}
