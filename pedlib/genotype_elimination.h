
inline int is_parent_child_genotype_compatible (int locus, int parent, int childSex,
						Genotype * pParentGeno,
						Genotype * pChildGeno);
inline int is_parent_child_allele_compatible (int alleleSetLen,
					      int *pParentAlleleSet,
					      int parent, int childSex,
					      Genotype * pChildGenotype);
