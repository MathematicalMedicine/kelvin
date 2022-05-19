/* Copyright (C) 2014, 2022 Mathematical Medicine LLC
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include "model_type.h"

double deNormalizeMean (ModelType *mt, double mean)
{
  if (mt->trait == DT || mt->distrib == QT_FUNCTION_CHI_SQUARE)
    return (mean);

  return (mean * mt->sd + mt->mean);
}

double deNormalizeStdev (ModelType *mt, double stdev)
{
  if (mt->trait == DT || mt->distrib == QT_FUNCTION_CHI_SQUARE)
    return (stdev);

  return (stdev * mt->sd);
}

