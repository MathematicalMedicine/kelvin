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
unsigned long cL[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 },    ///< Actual # calls to each instance of compute_likelihood
    eCL[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};     ///< Est. final # calls to each instance of compute_likelihood

FILE *fpCond = NULL;    ///< Conditional LR for genetic counseling, global due to likelihood.c write!
FILE *fpHet = NULL;     ///< Average HET LR file (Bayes Ratio file) pointer
FILE *fpPPL = NULL;     ///< PPL output file pointer
FILE *fpMOD = NULL;     // MOD and maximizing model information
FILE *fpIR = NULL;      ///< Intermediate results, used to go to stderr, normally dkelvin-only
FILE *fpDK = NULL;      // DCHURE detail file
FILE *fpDry = NULL;     ///< Dry-run statistics output for sizing estimation   

LambdaCell *pLambdaCell;
