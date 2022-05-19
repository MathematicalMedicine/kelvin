/* Copyright (C) 2008, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
	/* This part is invariant. */
	int i;
	void *libraryHandle;
	char functionName[256], libraryName[256];
	double (*pFunction)();

	for (i = 0; i < dLFunctionCount; i++) {
	  if ((i % 100) == 0) {
	    if (i != 0)
	      dlclose (libraryHandle);
	    sprintf (libraryName, "./%s_%04d.so", baseFunctionName, i);
	    libraryHandle = dlopen (libraryName, RTLD_LAZY|RTLD_LOCAL);
	    if (libraryHandle == NULL) {
	      fprintf (stderr, "dlopen() error [%s] for library file %s\n", dlerror(), libraryName);
	      exit (EXIT_FAILURE);
	    }
	  }
	  sprintf (functionName, "%s_%04d", baseFunctionName, i);
	  if ((pFunction = dlsym (libraryHandle, functionName)) == NULL) {
	    fprintf (stderr, "dlsym() error [%s] for function %s in library file %s\n", dlerror(),
		     functionName, libraryName);
	    exit (EXIT_FAILURE);
	  }
	  pFunction(V, S, P, F);
	}
	dlclose (libraryHandle);
