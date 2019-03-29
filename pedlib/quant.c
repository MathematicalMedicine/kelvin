
v_double
compute_quant (LOCI * loc, v_double * phenotype, int lindex, int right_index)
{
  v_double v_1;
  v_double penetrance;
  v_double v_penetrance;
  int j, k, num_traits;
  v_boolean unknown_quant;
  V_PAP = V_FALSE;
  penetrance = 0.0;
  num_traits = loc->num_liab_class;
  unknown_quant = V_TRUE;
  for (j = 0; j < num_traits; j++)
    {
      if (phenotype[j] > v_unknown_quant || phenotype[j] < v_unknown_quant)
	unknown_quant = V_FALSE;
    }
  if (unknown_quant == V_TRUE)
    {
      v_num_unknown++;
      return (1.0);
    }

  for (j = 0; j < num_traits; j++)
    {
#ifdef V_DEBUG
      if (v_debug_genotype['3'] == V_TRUE)
	fprintf (v_debugfile, cnvrt_format ("\n phen->x[%d] = %f"), j,
		 phenotype[j]);
#endif
      for (k = 0; k < num_traits; k++)
	{
	  if (lindex == right_index)
	    {
	      /*
	         v_1 = (phenotype[j] - loc->q_means[j][lindex][right_index]);
	         v_1 *=(phenotype[k] - loc->q_means[k][lindex][right_index]);
	       */
	      v_1 = (phenotype[j] - loc->q_means[j][lindex + right_index]);
#ifdef V_DEBUG
	      if (v_debug_genotype['3'] == V_TRUE)
		fprintf (v_debugfile, cnvrt_format ("\n v1 = %f"), v_1);
#endif
	      v_1 *= (phenotype[k] - loc->q_means[k][lindex + right_index]);
#ifdef V_DEBUG
	      if (v_debug_genotype['3'] == V_TRUE)
		fprintf (v_debugfile, cnvrt_format ("\n v1*v1 = %f"), v_1);
#endif
	      penetrance += loc->q_variance[j][k] * v_1;
#ifdef V_DEBUG
	      if (v_debug_genotype['3'] == V_TRUE)
		fprintf (v_debugfile, cnvrt_format ("\n v1*v1/sigma^2 = %f"),
			 penetrance);
#endif
	      if (V_PAP == V_TRUE)
		{
		  penetrance *= loc->q_variance[j][k];
		}
	      /*
	         penetrance +=loc->q_variance[j][k]*
	         (phenotype[j] - loc->q_means[j][lindex][right_index])*(phenotype[k] - loc->q_means[k][lindex][right_index]);
	       */
	    }
	  else
	    {
	      /*
	         v_1 = (phenotype[j] - loc->q_means[j][lindex][right_index]);
	         v_1 *=(phenotype[k] - loc->q_means[k][lindex][right_index]);
	       */
	      v_1 = (phenotype[j] - loc->q_means[j][lindex + right_index]);
#ifdef V_DEBUG
	      if (v_debug_genotype['3'] == V_TRUE)
		fprintf (v_debugfile, cnvrt_format ("\n v1 = %f"), v_1);
#endif
	      v_1 *= (phenotype[k] - loc->q_means[k][lindex + right_index]);
#ifdef V_DEBUG
	      if (v_debug_genotype['3'] == V_TRUE)
		fprintf (v_debugfile, cnvrt_format ("\n v1*v1 = %f"), v_1);
#endif
	      penetrance += loc->het_multiplier * loc->q_variance[j][k] * v_1;
#ifdef V_DEBUG
	      if (v_debug_genotype['3'] == V_TRUE)
		fprintf (v_debugfile, cnvrt_format ("\n v1*v1/sigma^2 = %f"),
			 penetrance);
#endif
	      if (V_PAP == V_TRUE)
		{
		  penetrance *= loc->q_variance[j][k];
		}
	      /* 
	         penetrance +=loc->het_multiplier*loc->q_variance[j][k]*
	         (phenotype[j] - loc->q_means[j][lindex][right_index])*(phenotype[k] - loc->q_means[k][lindex][right_index]);
	       */
	    }
#ifdef V_DEBUG
	  if (v_debug_genotype['3'] == V_TRUE)
	    {
	      fprintf (v_debugfile, cnvrt_format ("\n q_means[%d] = %f"), j,
		       loc->q_means[k][lindex + right_index]);
	      fprintf (v_debugfile,
		       cnvrt_format ("\n (variance )[%d][%d] = %10.8f"), j, k,
		       loc->q_variance[j][k]);
	    }
#endif
	}

    }

#ifdef V_DEBUG
  if (v_debug_genotype['3'] == V_TRUE)
    fprintf (v_debugfile, cnvrt_format ("\n  determ = %10.8f"), covar_determ);
#endif

  v_penetrance = -0.5 * penetrance;
  penetrance = covar_determ * dexp (v_penetrance);
  if (V_PAP == V_TRUE)
    {
      penetrance *= covar_determ;
    }


#ifdef V_DEBUG
  if (v_debug_genotype['3'] == V_TRUE)
    {
      fprintf (v_debugfile, cnvrt_format ("    exponent = %10.8f"),
	       v_penetrance);
      fprintf (v_debugfile, cnvrt_format ("\n  dexp(pen) = %10.8f"),
	       dexp (v_penetrance));
      fprintf (v_debugfile, cnvrt_format ("\n 1/mult = %10.8f"),
	       1 / (covar_determ * covar_determ));
      fprintf (v_debugfile, cnvrt_format ("\n  Val = %10.8f"), penetrance);
    }
#endif

  if (lindex != right_index)
    penetrance *= loc->cov_adj;

  if (V_PAP == V_TRUE)
    {
      penetrance = penetrance / (dsqrt (2.0 * V_PI));
    }
#ifdef V_DEBUG
  if (v_debug_genotype['3'] == V_TRUE)
    fprintf (v_debugfile, cnvrt_format ("\n PEN Val = %10.8f\n"), penetrance);
#endif

  return (penetrance);
}
