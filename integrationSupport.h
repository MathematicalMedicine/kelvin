/* Copyright (C) 2008, 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#define SCALE_RESERVE 50
#define checkpt() fprintf (stderr, "Checkpoint at line %d of file \"%s\"\n",__LINE__,__FILE__)
void compute_hlod_2p_dt (double x[], double *f, int *scale);
int kelvin_dcuhre_integrate (double *integral, double *abserr, double, int *);
void compute_hlod_mp_dt (double x[], double *f, int *scale);
void compute_hlod_2p_qt (double x[], double *f, int *scale);
void compute_hlod_mp_qt (double x[], double *f, int *scale);
void integrateMain();



