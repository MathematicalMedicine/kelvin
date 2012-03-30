#define SCALE_RESERVE 50
#define checkpt() fprintf (stderr, "Checkpoint at line %d of file \"%s\"\n",__LINE__,__FILE__)
void compute_hlod_2p_dt (double x[], double *f, int *scale, double []);
int kelvin_dcuhre_integrate (double *integral, double *abserr, double, int *);
void compute_hlod_mp_dt (double x[], double *f, int *scale,double []);
void compute_hlod_2p_qt (double x[], double *f, int *scale,double []);
void compute_hlod_mp_qt (double x[], double *f, int *scale,double []);
void integrateMain();



