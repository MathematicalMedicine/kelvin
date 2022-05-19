/* Copyright (C) 2006, 2010, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include "model.h"

void initializeDefaults ();
void readConfigFile (char *config);
void parseCommandLine (int argc, char *argv[]);
void validateConfig ();
void fillConfigDefaults (ModelRange *modelRange, ModelOptions *modelOptions, ModelType *modelType);
void finishConfig (ModelRange *modelRange, ModelType *modelType);
