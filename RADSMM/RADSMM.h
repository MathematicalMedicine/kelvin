 /* Random Access Data Storage for MLip (Multiple Model) (RADSMM) */
  
#define MAX_DATA_FILES 20
#define RADSMM_descr_length 256
        
/* Some constants, usually only used for range checking */
#define RADSMM_MAX_markers 200
#define RADSMM_MAX_pedigrees 1000
#define RADSMM_MAX_thetas  84000
#define RADSMM_MAX_penetrances 176750 
#define RADSMM_MAX_qmodels 80000 
#define RADSMM_MAX_diseqs 80000 
#define RADSMM_MAX_geneFreqs 100
#define RADSMM_MAX_markerlabel_length 48
#define RADSMM_MAX_pedigreelabel_length 80
#define RADSMM_MAX_liability_classes 12
#define RADSMM_MAX_diseq_params 400

typedef int    RADSMM_pedigree_type;
typedef float  RADSMM_marker_type;
/* The following three were originally double */
typedef float  RADSMM_theta_type;
typedef float  RADSMM_geneFreq_type;
typedef float  RADSMM_diseq_type;
typedef struct  penet {
        float   penetrance[3];  
} RADSMM_penetrance_type;
typedef struct  qtype {
        float   means[3];  
        float   variances[3];  
} RADSMM_quantitative_type;

typedef struct RADSMM_proto_header {   
        long             marker_count;
        RADSMM_marker_type	*marker_list;
        long             pedigree_count;
        RADSMM_pedigree_type	*pedigree_list;
        long             theta_count;
        RADSMM_theta_type	*theta_list;
        long             LC_count;
        long             penetrance_count;
        RADSMM_penetrance_type	*penetrance_list[RADSMM_MAX_liability_classes];
        long             qmodel_count;
        RADSMM_quantitative_type	*qmodel_list;
        long             diseq_count;
	RADSMM_diseq_type 	*diseq_list;
        long             geneFreq_count;
        RADSMM_geneFreq_type	*geneFreq_list;
	long		markerlabel_size;  /* max size of a label string */
	char		*markerlabel_list;
	long		pedigreelabel_size;
	char		*pedigreelabel_list;
/* Misc Stuff */
	int		version;
	int		subversion;
        char            element_data_type;  /* 'D'=double, 'F'=float */
        char            theta_matrix_type;  /* 'G'=grid, 'D'=diagonal */
        char    	description[RADSMM_descr_length];
        char    	date_string[17];
	char		ordering;   /* order of parameter presidence */
        char            marker_type;  /* '2' for 2 point, 'M' for multi-point */
	char		model_type;  /* 'D'icot or 'Q'uant */
        char		use_Diseq;  /* Use lambda */
        int             number_of_files;  /* really the number of extra files */
/* Internal use */
	int	checksum_type;	  /* not used yet.. */
        long	start_of_data;	  /* Where the header ends and data begins in file */
        long    pedigree_offset;  /* these next 3 vars are not really used */
        long    marker_offset;
        long    theta_offset;
        long    penetrance_offset;
        long    geneFreq_offset;
        long    qmodel_offset;
        long    diseq_offset;
	long	chunk_size;	  /* 8 if double, 4 if floats */
	long	chunks_per_file;   /* chunks per file (not including header) */
/* not part of the file */
	int	debug_level;		/* bitmapped, see doc */
        char	open_flag;    	  /* 'O' = open and ready to read/write */
        int	fp[MAX_DATA_FILES];

} RADSMM_header_type;

/* Non-numerical data values */
#define RADSMM_EMPTY             3.29E+38	/* Initialized state */
#define RADSMM_IGNORED           3.28E+38	/* Purposefully left blank */
        
#define RADSMM_Infinity          3.40E+38
#define RADSMM_Negative_Infinity 3.39E38
#define RADSMM_Not_Possible      3.38E38
#define RADSMM_Not_a_Number      3.37E38
#define RADSMM_not_data_limit    3.20E38 

/* error Codes */
#define RADSMM_ERROR_success 0
#define RADSMM_ERROR_badindex -1
#define RADSMM_ERROR_badparam -2
#define RADSMM_ERROR_badpointer -3
#define RADSMM_ERROR_lseek -4
#define RADSMM_ERROR_fileopen -5 
#define RADSMM_ERROR_writing -6
#define RADSMM_ERROR_locking -7
#define RADSMM_ERROR_malloc -8
#define RADSMM_ERROR_reading -9
#define RADSMM_ERROR_internal -10
#define RADSMM_ERROR_not_open -11
#define RADSMM_ERROR_already_open -12
#define RADSMM_ERROR_file_header -13
#define RADSMM_ERROR_value_not_in_list -15
#define RADSMM_ERROR_writeover_valid_data -16
#define RADSMM_ERROR_outofrange -17
#define RADSMM_ERROR_wrong_model -18

/* lets define the subroutines */
extern int RADSMM_setup_init( RADSMM_header_type *header, int debug_level );
extern int RADSMM_setup_marker( RADSMM_header_type *header, RADSMM_marker_type marker_list[], long marker_count );
extern int RADSMM_setup_pedigree( RADSMM_header_type *header, int pedigree_list[], long pedigree_count );
extern int RADSMM_setup_theta( RADSMM_header_type *header, double theta_list[], long theta_count, char theta_matrix_type );
extern int RADSMM_setup_LC( RADSMM_header_type *header, long count );
extern int RADSMM_setup_penetrance( RADSMM_header_type *header, int LC_ndx, float plist1[],float plist2[],float plist3[], long penetrance_count );
extern int RADSMM_setup_geneFreq( RADSMM_header_type *header, double genefreq[], long genefreq_count );
extern int RADSMM_setup_qmodel( RADSMM_header_type *header, RADSMM_quantitative_type list[], long count );
extern int RADSMM_setup_diseq( RADSMM_header_type *header, RADSMM_diseq_type lambda_list[], long lambda_count );
extern int RADSMM_setup_data( RADSMM_header_type *header, char Data_type, int checksum_type );
extern int RADSMM_setup_comments( RADSMM_header_type *header, char *comment );
extern int RADSMM_setup_ordering( RADSMM_header_type *header, char ordering );
extern int RADSMM_setup_type( RADSMM_header_type *header, char Point_Type, char model_type, char use_diseq );
extern int RADSMM_file_size( RADSMM_header_type *header, double *filesize, int *number_of_files );
extern int RADSMM_create_file( RADSMM_header_type *header,  char *filename, int mode );
extern int RADSMM_open_file( RADSMM_header_type *header,  char *filename, char open_mode );
extern int RADSMM_close_file( RADSMM_header_type *header );
extern long RADSMM_index_marker( RADSMM_header_type *header, double location );
extern long RADSMM_index_pedigree( RADSMM_header_type *header, int pedid );
extern long RADSMM_index_theta( RADSMM_header_type *header, double male_theta, double female_theta );
extern long RADSMM_index_penetrance( RADSMM_header_type *header, int LC, float penetrance1, float penetrance2, float penetrance3 );
extern long RADSMM_index_genefreq( RADSMM_header_type *header, double genefreq );
extern long RADSMM_index_diseq( RADSMM_header_type *header, double lambda );

extern int RADSMM_read_data( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
	long genefreq_index, long penet_index, long qmodel_index, long diseq_index, double *data );

extern int RADSMM_write_data( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
	long genefreq_index, long penet_index, long qmodel_index, long diseq_index, double *data );

extern int RADSMM_read_float( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
	long genefreq_index, long penet_index, long qmodel_index, long diseq_index, float *data );

extern int RADSMM_write_float( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
	long genefreq_index, long penet_index, long qmodel_index, long diseq_index, float *data );

extern int RADSMM_read_list_float( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
	long genefreq_index, long penet_index, long qmodel_index, long diseq_index, float data[], long cnt );

extern int RADSMM_write_list_float( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
	long genefreq_index, long penet_index, long qmodel_index, long diseq_index, float data[], long cnt );

extern int RADSMM_read_list_double(RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
	long genefreq_index, long penet_index, long qmodel_index, long diseq_index, double data[], long cnt );


extern int RADSMM_write_list_double( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index, 
	long genefreq_index, long penet_index, long qmodel_index, long diseq_index, double data[], long cnt );

extern int RADSMM_sync( RADSMM_header_type *header );
extern int RADSMM_setup_markerlabel( RADSMM_header_type *header, int label_length );
extern int RADSMM_set_markerlabel( RADSMM_header_type *header, int markerindex,
                            char *label );
extern int RADSMM_setup_pedigreelabel( RADSMM_header_type *header, int label_length );
extern int RADSMM_set_pedigreelabel( RADSMM_header_type *header, int pedigree_index,
                            char *label );
extern char  *RADSMM_get_markerlabel( RADSMM_header_type *header, int marker_index ); 
extern char  *RADSMM_get_pedigreelabel( RADSMM_header_type *header,
                                 int pedigree_index );
extern long RADSMM_index_markerlabel( RADSMM_header_type *header, char *marker );
extern long RADSMM_index_pedigreelabel( RADSMM_header_type *header,
                                 char *pedigreelabel );

/* internal subroutine */
extern int RADSMM_seek( RADSMM_header_type *header, long ped_index, long marker_index, long theta_index,
   long genefreq_index, long penet_index, long qmodel_index, long diseq_index, int *fpndx );

extern int RADSMM_range_check(   RADSMM_header_type *header, long ped_index, long marker_index, long theta_index,
   long genefreq_index, long penet_index, long qmodel_index, long diseq_index, long ndxcnt );

extern int RADSMM_open_others( char filename[], int number, int flag, int open_mode  );

