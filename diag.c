#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>
#include "sw.h"

int main (int argc, char *argv[]) {
  int i, segid;
  int *envDiagLevel;

  segid = atoi (argv[1]);
  INFO ("Attaching to segment ID %d", segid);
  if ((envDiagLevel = shmat (segid, NULL, 0)) == ((void *) -1))
    ERROR ("Cannot attach shared memory segment for diagnostics, %s", strerror(errno));
  for (i=0; i<MAX_DIAG_FACILITY; i++) {
    fprintf (stderr, "envDiagLevel[%d] was %d, now ", i, envDiagLevel[i]);
    envDiagLevel[i] = 0;
    fprintf (stderr, "%d\n", envDiagLevel[i]);
  }
  exit (EXIT_SUCCESS);
}
