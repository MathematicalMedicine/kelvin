
/**********************************************************************
 * Copyright 2010, Nationwide Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

// NOTE: functions are no longer inlined due to a conflict between clang and
// gcc's interpretations of the "inline" keyword.
/*inline */int is_parent_child_genotype_compatible (int locus, int parent,
						int childSex,
						Genotype * pParentGeno,
						Genotype * pChildGeno);
/*inline */int is_parent_child_allele_compatible (int alleleSetLen,
					      int *pParentAlleleSet,
					      int parent, int childSex,
					      Genotype * pChildGenotype);
/*inline */int isHet (Genotype * pGeno);
