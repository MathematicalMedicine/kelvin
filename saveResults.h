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

void dump_lDT (double **);
int saveTrait (int, char *, double **);
int restoreTrait (int, char *, double **);
int saveMarker (char *, int, int, char **, double *);
int restoreMarker (char *, int, int, char **, double *);
int saveAlternative (char *, int, double, double **);
int restoreAlternative (char *, int, double, double **);
