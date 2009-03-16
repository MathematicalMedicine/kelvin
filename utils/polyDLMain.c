
int main (int argc, char *argv[]) {
  int i;
  struct polynomial **V;
  if ((argc - 1) != baseFunctionArgs) {
    fprintf (stderr, "Function %s requires %d floating arguments; %d provided on command line\n",
	     baseFunctionName, baseFunctionArgs, argc - 1);
    exit (EXIT_FAILURE);
  }
  V = (struct polynomial **) malloc ((sizeof (struct polynomial *)) * baseFunctionArgs);
  printf ("Calling %s with ", baseFunctionName);
  for (i=1; i<argc; i++) {
    if (i != 1) printf (", ");
    V[i - 1] = (struct polynomial *) malloc (sizeof (struct polynomial));
    V[i - 1]->value = atof (argv[i]);
    printf ("%g", V[i - 1]->value);
  }
  printf ("\n");
  printf ("= %g\n", (baseFunction)(1, V));
  exit (EXIT_SUCCESS);
}
