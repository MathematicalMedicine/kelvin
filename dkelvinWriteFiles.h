/* Copyright (C) 2009, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
void dk_write2ptBRHeader (int loc1, int loc2);
void dk_write2ptBRData (double dprimevalue, double theta1,double theta2,double integral, int max_scale);
void dk_writeMPBRHeader ();
void dk_writeMPBRData (int posIdx, float traitPos, double ppl, double br, int max_scale);
void dk_writeMPMODHeader ();
void dk_writeMPMODData (int posIdx, float traitPos, double value, st_DKMaxModel *model);
void dk_write2ptMODHeader ();
void dk_write2ptMODData (char *description, double value, st_DKMaxModel *model);
void dk_copyMaxModel (double *arr, st_DKMaxModel *max, int num);
void dk_copyMaxModel2 (st_DKMaxModel *dest, st_DKMaxModel *cur);
