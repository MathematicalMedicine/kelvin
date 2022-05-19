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
void writePPLFileHeader ();
void writePPLFileDetail (int dprime0Idx);
void write2ptBRFile (int loc1, int loc2);
void writeMPBRFileHeader ();
void writeMPBRFileDetail (int posIdx, float traitPos, float ppl, double avgLR);
void writeMPMODFileHeader ();
void writeMPMODFileDetail (int posIdx, float traitPos);
void writeMaximizingModel (char *modelDescription, double myMOD, int myDPrimeIdx,
			   int myThetaIdx);
void write2ptMODFile (int loc1, int loc2, int dprime0Idx);
void writeSurfaceFileHeader ();






