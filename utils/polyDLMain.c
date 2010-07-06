#ifdef NOTKELVIN

int main () {
  int i;
  double vD;
  Polynomial **V;
  V = (Polynomial **) malloc ((sizeof (Polynomial *)) * baseFunctionArgs);
  for (i=0; i<baseFunctionArgs; i++) {
    V[i] = (Polynomial *) malloc (sizeof (Polynomial));
    V[i]->e.v = (struct variablePoly *) malloc (sizeof (struct variablePoly));
    printf ("Name for variable number %d (no spaces or punctuation): ", i);
    scanf ("%s", V[i]->e.v->vName);
    printf ("Floating value for variable named %s (floating number only!): ", V[i]->e.v->vName);
    scanf ("%f", &V[i]->value);
  }
  printf ("\n");
  printf ("= %g\n", (baseFunction)(1, V));
  exit (EXIT_SUCCESS);
}

#else

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

#endif
