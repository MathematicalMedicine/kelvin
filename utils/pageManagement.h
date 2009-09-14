#ifndef __PAGE_MANAGEMENT__
#define __PAGE_MANAGEMENT__

/**
@file pageManagement.h

  @see pageManagement.c
  
  Copyright &copy; 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$

*/

void setupSegvHandler (void);
void *allocatePages (int objectSizeInBytes);
void allowReadOnly (void *object, int objectSizeInBytes);
void allowReadWrite (void *object, int objectSizeInBytes);

#endif
