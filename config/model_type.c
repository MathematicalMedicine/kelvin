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

