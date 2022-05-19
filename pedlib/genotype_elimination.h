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
