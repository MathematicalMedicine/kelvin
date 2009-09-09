#include <unistd.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>

void *allocatePages (int objectSizeInBytes)
{
  void *pageStart;
  int pageSize, pageCount;

  pageSize = getpagesize();
  pageCount = objectSizeInBytes / pageSize + 1;

  printf ("There are %d pages from %d bytes and %d bytes/page\n", pageCount,
	  objectSizeInBytes, pageSize);

  /* Allocate discrete pages of accessable memory for this structure. */
  if (NULL == (pageStart =
	       mmap(0 /* Hinted start */,
		    (pageCount * pageSize) /* Size in bytes */, 
		    PROT_READ|PROT_WRITE /* Protection */, 
		    MAP_ANON|MAP_PRIVATE /* Flags */,
		    -1 /* FD, ignored given MAP_ANON */,
		    0 /* Offset, ignored given MAP_ANON */))) {
    fprintf (stderr, "Memory page allocation for %d-byte structure failed!\n",
	     objectSizeInBytes);
    exit (EXIT_FAILURE);
  }
  return pageStart;
}

void allowReadOnly (void *pageStart, int objectSizeInBytes)
{
  int pageSize, pageCount;

  pageSize = getpagesize();
  pageCount = objectSizeInBytes / pageSize + 1;

  if (mprotect(pageStart, pageCount, PROT_READ) == -1)
    perror("mprotect");
}

void allowReadWrite (void *pageStart, int objectSizeInBytes)
{
  int pageSize, pageCount;

  pageSize = getpagesize();
  pageCount = objectSizeInBytes / pageSize + 1;

  if (mprotect(pageStart, pageCount, PROT_READ|PROT_WRITE) == -1)
    perror("mprotect");
}

#ifdef MAIN

//#include "utils/pageManagement.h"

int main (int argc, char * argv[])
{
  struct little_one {
    int beginning;
    int ending;
  };
  struct big_one {
    int beginning;
    void *middle[1024];
    int ending;
  };
  struct little_one *littleOne;
  struct big_one *bigOne;

  littleOne = (struct little_one *) allocatePages (sizeof (struct little_one));
  bigOne = (struct big_one *) allocatePages (sizeof (struct big_one));

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

  fprintf (stdout, "We should now get a segmentation fault attempting to write to the protected structure:\n");

  allowReadOnly (littleOne, sizeof (struct little_one));

  littleOne->beginning = 3;
}

#endif
