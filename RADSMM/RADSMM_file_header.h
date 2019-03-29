/* this define the actual file header */
/* RADSMM Version 1.0 */
typedef struct RADSMM_fileheader {   
	char		cookie[4];
	int		version;
	int		subversion;		/* subversion means compatible change */
	long		start_of_data;		/* Offset of header & indexs */

	long		marker_count;
        long            marker_offset;

	long		pedigree_count;
	long		pedigree_offset;

	long		theta_count;
	long		theta_offset;
	char		theta_matrix_type;  /* 'G'=grid, 'D'=diagonal */

	char		padding1[3];

	long		penetrance_count;
        long            LC_count;
	long		penetrance_offset;

	long		qmodel_count;   /* Q model */
	long		qmodel_offset;

	long		diseq_count;
	long		diseq_offset;

	long		geneFreq_count;
	long		geneFreq_offset;

	long		markerlabel_size;
	long		markerlabel_offset;

	long		pedigreelabel_size;
	long		pedigreelabel_offset;

	char		element_data_type;	/* 'F' or 'D' */
	char		model_type; /* 'D'icot or 'Q'uant */
	char		marker_type; /* '2' for 2 point or 'M' for Multipoint */
        char 		use_Diseq; /* 'N' for no or anything else to enable */

	long		chunks_per_file;	/* for multiple files */
	int		number_of_files;
	char		ordering;
	char		padding2[32];		/* add new vars here. */

        char    	date_string[17];
        char    	description[64];
} RADSMM_file_header_type;
/* after this header the indexs and labels are stored */


/* a few constansts */

#define RADSMM_cookie "RDMM"   /* must be 4 characters */
