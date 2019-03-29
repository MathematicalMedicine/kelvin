/**********************************************************************
 * Copyright 2007, Columbus Children's Research Institute.  
 * All rights reserved.
 * Permission is hereby given to use this software 
 * for non-profit educational purposes only.
 **********************************************************************/

inline int is_parent_child_genotype_compatible (int locus, int parent, int childSex,
						Genotype * pParentGeno,
						Genotype * pChildGeno);
inline int is_parent_child_allele_compatible (int alleleSetLen,
					      int *pParentAlleleSet,
					      int parent, int childSex,
					      Genotype * pChildGenotype);
