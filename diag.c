/* Copyright (C) 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include <sys/ipc.h>
#include <sys/shm.h>
#include <errno.h>
#include "utils/sw.h"

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
