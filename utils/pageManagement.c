/**
@file pageManagement.c

  Page-at-a-time replacement for (m|c)alloc that allows us to read-lock global variables.

  Since we're not going to use C++ and therefore can't have restricted accessors for the
  global variables that would otherwise be objects, we do the next best thing by allocating
  space for these global variables a page and a time and using mprotect() to mark them
  as read-only once they're initialized.

  A good example of usage is in config/config.c, where the signal handler is setup and at
  least some configuration-related global variables allocated and read-locked. This is, of
  course, primarily a development tool to identify bad behavior, so it could be conditionalized
  away for production runs.

  Copyright &copy; 2009, Nationwide Children's Research Institute.  All
  rights reserved.  Permission is hereby given to use this software
  for non-profit educational purposes only.

  @version $Id$
*/
  
#include <unistd.h>
#include <sys/mman.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

/*
  The signals for violations of memory protection depend upon 
  OS lineage. Somehow I felt that this should be __BSD__, but
  nothing like that exists.
*/
#ifdef __APPLE__
  #define FAULT_NAME SIGBUS
#else
  #define FAULT_NAME SIGSEGV
#endif

static void segvHandler (int signum, siginfo_t *sigi, void *unused)
{
  if (sigi->si_signo == FAULT_NAME && sigi->si_code == SEGV_ACCERR) {
    fprintf (stderr, "Protection error referencing stucture at address %lx\n",
	     (unsigned long) sigi->si_addr);
    abort ();
  }
}


/**

  Setup the handler for bus errors for segmentation violations
  that we might cause by read-protecting our globals.

  @author Bill Valentine-Cooper - overall content.
  @par Reviewers:

  @par Global Inputs

  none

  @par Global Outputs

  none

  @return Nothing

*/
void setupSegvHandler (void)
{
  struct sigaction segvAction;
  segvAction.sa_handler = (void *) segvHandler;
  segvAction.sa_flags = SA_RESTART | SA_SIGINFO | SA_NODEFER | SA_RESETHAND;
  sigemptyset (&segvAction.sa_mask);
  if (sigaction (FAULT_NAME, &segvAction, NULL)) {
    perror ("sigaction");
    exit (EXIT_FAILURE);
  }
}

/**

  Allocate a series of contiguous pages large enough to hold the number
  of bytes specified and return a pointer to the start address.

  @author Bill Valentine-Cooper - overall content.

  @par Global Inputs

  Process address space

  @par Global Outputs

  A little less free address space

  @return a pointer to the start of the pages allocated, exits with EXIT_FAILURE if the allocation fails.

*/
void *allocatePages (
		     int objectSizeInBytes ///< Number of bytes to allocate
		     )
{
  void *pageStart;
  int pageSize, pageCount;

  pageSize = getpagesize ();
  pageCount = objectSizeInBytes / pageSize + 1;

  /* Allocate discrete pages of accessable memory for this structure. */
  if (NULL == (pageStart = mmap (0 /* Hinted start */ ,
              (pageCount * pageSize) /* Size in bytes */ ,
              PROT_READ | PROT_WRITE /* Protection */ ,
              MAP_ANON | MAP_PRIVATE /* Flags */ ,
              -1 /* FD, ignored given MAP_ANON */ ,
              0 /* Offset, ignored given MAP_ANON */ ))) {
    fprintf (stderr, "Memory page allocation for %d-byte structure failed!\n", objectSizeInBytes);
    exit (EXIT_FAILURE);
  }
  return pageStart;
}

/**

  Change protection for the page(es) specified to allow only read references.

  @author Bill Valentine-Cooper - overall content.

  @par Global Inputs

  None

  @par Global Outputs

  None

  @return none, exits with EXIT_FAILURE if the mprotect fails.

*/
void allowReadOnly (
		    void *pageStart, ///< Starting address of pages to be protected
		    int objectSizeInBytes ///< Size of region to protect in bytes, will convert to pages
		    )
{
  int pageSize, pageCount;

  pageSize = getpagesize ();
  pageCount = objectSizeInBytes / pageSize + 1;

  if (mprotect (pageStart, pageCount, PROT_READ) == -1) {
    perror ("mprotect");
    exit (EXIT_FAILURE);
  }
}

/**

  Change protection for the page(es) specified to allow read and write access.

  @author Bill Valentine-Cooper - overall content.

  @par Global Inputs

  None

  @par Global Outputs

  None

  @return none, exits with EXIT_FAILURE if the mprotect fails.

*/
void allowReadWrite (
		    void *pageStart, ///< Starting address of pages to be protected
		    int objectSizeInBytes ///< Size of region to protect in bytes, will convert to pages
		    )
{
  int pageSize, pageCount;

  pageSize = getpagesize ();
  pageCount = objectSizeInBytes / pageSize + 1;

  if (mprotect (pageStart, pageCount, PROT_READ | PROT_WRITE) == -1)
    perror ("mprotect");
}

#ifdef MAIN

//#include "utils/pageManagement.h"

int main (int argc, char *argv[])
{
  struct little_one
  {
    int beginning;
    int ending;
  };
  struct big_one
  {
    int beginning;
    void *middle[1024];
    int ending;
  };
  struct little_one *littleOne;
  struct big_one *bigOne;

  setupSegvHandler ();

  littleOne = (struct little_one *) allocatePages (sizeof (struct little_one));
  fprintf (stderr, "littleOne starts at %lx\n", littleOne);
  bigOne = (struct big_one *) allocatePages (sizeof (struct big_one));
  fprintf (stderr, "bigOne starts at %lx\n", bigOne);

  littleOne->beginning = 2.0;
  littleOne->ending = 3.0;

  bigOne->beginning = 12.0;
  bigOne->ending = 13.0;

  allowReadOnly (littleOne, sizeof (struct little_one));
  allowReadOnly (bigOne, sizeof (struct big_one));

  fprintf (stdout, "littleOne->beginning is %d\n", littleOne->beginning);
  fprintf (stdout, "littleOne->ending is %d\n", littleOne->ending);
  fprintf (stdout, "bigOne->beginning is %d\n", bigOne->beginning);
  fprintf (stdout, "bigOne->ending is %d\n", bigOne->ending);

  allowReadWrite (littleOne, sizeof (struct little_one));
  littleOne->beginning = 2;

  fprintf (stdout, "littleOne->beginning is %d\n", littleOne->beginning);

  fprintf (stdout, "We should now get a 'segmentation fault'/'bus error' attempting to write to the protected structure:\n");

  allowReadOnly (littleOne, sizeof (struct little_one));

  littleOne->ending = 3;
}

#endif
