#include <unistd.h>
#include <sys/mman.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>


/* Now for the grotty machine-dependent stuff.  */
#ifdef __APPLE__
#define SIGSEGV SIGBUS
#endif

#if defined (__linux__) && (defined (__powerpc__) || defined (__i386__))

  /* For a PowerPC or i386 box running Linux. */

  #if defined (__powerpc__)
    #define EXTRA_ARGS , struct sigcontext_struct *sigc
    #define FAULT_ADDRESS ((void *)sigc->regs->dar)
    #define SIGNAL_OK 					\
      (sigc->signal == SIGSEGV 				\
       && sigc->regs->trap == 0x00300			\
       && (sigc->regs->dsisr & 0x02000000) != 0)
  #endif
  #if defined (__i386__)
    #define EXTRA_ARGS , struct sigcontext sigc
    #define FAULT_ADDRESS ((void*)(sigc.cr2))
    #define SIGNAL_OK (sigc.trapno == 14)
  #endif

  static struct sigaction segv_action;
  #define CONTINUE sigaction (SIGSEGV, &segv_action, NULL); return

  #define SETUP_HANDLER(handler)						\
    do {									\
      sigaction (SIGSEGV, NULL, &segv_action);				\
      segv_action.sa_handler = (sig_t)(handler);				\
      /* We want:								\
         system calls to resume if interrupted;				\
         segfaults inside the handler to be trapped, by the default	\
           technique (that is, a core dump).  */				\
      segv_action.sa_flags =						\
         SA_RESTART | SA_NODEFER | SA_RESETHAND;				\
      sigaction (SIGSEGV, &segv_action, NULL);				\
    } while (0)

#else

  /* For everything but a PowerPC or i386 box running Linux. */

  #define EXTRA_ARGS , siginfo_t *sigi, void *unused
  #define SIGNAL_OK (sigi->si_signo == SIGSEGV && sigi->si_code == SEGV_ACCERR)
  #define FAULT_ADDRESS (sigi->si_addr)

  static struct sigaction segv_action;
  #define CONTINUE sigaction (SIGSEGV, &segv_action, NULL); return

  #define SETUP_HANDLER(handler)						\
    do {									\
      sigaction (SIGSEGV, NULL, &segv_action);				\
      segv_action.sa_sigaction = (handler);				\
      /* We want:								\
         system calls to resume if interrupted;				\
         to use the three-argument version of the signal handler; and	\
         segfaults inside the handler to be trapped, by the default	\
           technique (that is, a core dump).  */				\
      segv_action.sa_flags =						\
         SA_RESTART | SA_SIGINFO | SA_NODEFER | SA_RESETHAND;		\
      sigaction (SIGSEGV, &segv_action, NULL);				\
    } while (0)

#endif

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

static void
segvHandler (int signum EXTRA_ARGS)
{
  fprintf (stderr, "Woot!\n");
  //  if (! SIGNAL_OK
  //      || ! maybe_mark_address (FAULT_ADDRESS))
  //    {
      /* For ease of debugging.  */
      static volatile void *fa;
      fa = FAULT_ADDRESS;
      fprintf (stderr, "Fault address is %lx\n", fa);
      abort ();
      //    }
  CONTINUE;
}

static void
setup_segvHandler (void)
{
  //  page_size_g = (size_t) sysconf (_SC_PAGESIZE);
  //  lock = 0;
  SETUP_HANDLER (segvHandler);
}

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

  setup_segvHandler ();

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
