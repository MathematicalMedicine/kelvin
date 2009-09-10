#include <unistd.h>
#include <sys/mman.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
  #define FAULT_NAME SIGBUS
#else
  #define FAULT_NAME SIGSEGV
#endif

static void segvHandler (int signum, siginfo_t *sigi, void *unused)
{
  if (sigi->si_signo == FAULT_NAME && sigi->si_code == SEGV_ACCERR) {
    fprintf (stderr, "Protection error referencing stucture at address %lx\n",
	     sigi->si_addr);
    abort ();
  }
}

void setupSegvHandler (void)
{
  struct sigaction segvAction;
  segvAction.sa_handler = (void *) segvHandler;
  segvAction.sa_flags = SA_RESTART | SA_SIGINFO | SA_NODEFER | SA_RESETHAND;
  if (sigaction (FAULT_NAME, &segvAction, NULL)) {
    perror ("sigaction");
    exit (EXIT_FAILURE);
  }
}

void *allocatePages (int objectSizeInBytes)
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

void allowReadOnly (void *pageStart, int objectSizeInBytes)
{
  int pageSize, pageCount;

  pageSize = getpagesize ();
  pageCount = objectSizeInBytes / pageSize + 1;

  if (mprotect (pageStart, pageCount, PROT_READ) == -1)
    perror ("mprotect");
}

void allowReadWrite (void *pageStart, int objectSizeInBytes)
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
