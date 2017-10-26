#include "common.h"
#include <time.h>

#include "fftsaxs.h"
#include "index.h"

#include "mol2/atom_group.h"
#include "mol2/pdb.h"
#include "mol2/prms.h"
#include "mol2/vector.h"

#ifdef _MPI_
#include <mpi.h>
#endif

void print_usage(char *exe_name)
{
  fprintf(stderr,
          "Usage: %s MAPPING_PRM ATOMPRM REC LIG "
          "Z_START Z_END L\n",
          basename(exe_name));
  
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	char* map_path = argv[1];
	char* prm_path = argv[2];
	char* rec_path = argv[3];
	char* lig_path = argv[4];
	char* exp_path = argv[5];
	double z_beg   = atof(argv[6]);
	double z_end   = atof(argv[7]);
	int    L       = atoi(argv[8]);
	
	double z_step = 1.0;
	
	if (argc != 9) { print_usage("puresaxs"); }
	
	int myrank = 0;
	int nprocs = 1;
	double starttime, endtime;
	
#ifdef _MPI_
	if ((MPI_Init(&argc, &argv)) != MPI_SUCCESS)
	{
		(void) fprintf (stderr, "MPI_Init failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
	if (MPI_Comm_rank(MPI_COMM_WORLD, &myrank) != MPI_SUCCESS)
		MPI_Abort (MPI_COMM_WORLD, 1);
	if (MPI_Comm_size(MPI_COMM_WORLD, &nprocs) != MPI_SUCCESS)
		MPI_Abort (MPI_COMM_WORLD, 1);
	
	starttime = MPI_Wtime();
#else
	starttime = (double)clock() / CLOCKS_PER_SEC;
#endif
	
	
	int nbeta = L + 1;
	int znum_total = (int)round((z_end - z_beg) / z_step) + 1;

	int qnum = QNUM;
	double* qvals = sxs_mkarray(0.0, QMAX, qnum); 

	if (myrank == 0) { SXS_PRINTF("Reading parameters ...\n"); }
	struct mol_prms *prms = mol_prms_read(prm_path);
	
	if (myrank == 0) { SXS_PRINTF("Reading form-factors ...\n"); }
	struct saxs_form_factor_table *ff_table = default_ff_table(map_path);
	
//=================================================================
//=========================== rec =================================
	if (myrank == 0) { SXS_PRINTF("Reading receptor ...\n"); }
	
	struct mol_atom_group *pa1 = mol_read_pdb(rec_path);
	mol_atom_group_add_prms(pa1, prms);
	
	struct mol_vector3 coe;
	center_of_extrema(&coe, pa1);
	MOL_VEC_MULT_SCALAR(coe, coe, -1.0);
	mol_atom_group_translate(pa1, &coe);

	struct sxs_spf_full* A = atom_grp2spf(pa1, ff_table, qvals, qnum, L, 1);

//=================================================================
//=========================== lig =================================
	if (myrank == 0) { SXS_PRINTF("Reading ligand ...\n"); }
	
	struct mol_atom_group *pa2 = mol_read_pdb(lig_path);
	mol_atom_group_add_prms(pa2, prms);
	
	struct mol_vector3 com;
	centroid(&com, pa2);
	MOL_VEC_MULT_SCALAR(com, com, -1.0);
	mol_atom_group_translate(pa2, &com);
	
	struct sxs_spf_full* B = atom_grp2spf(pa2, ff_table, qvals, qnum, L, 1);
	
//=================================================================
//=================== experiment ==================================
	if (myrank == 0) { SXS_PRINTF("Reading experiment ...\n"); }
	
	struct sxs_profile* exp_profile = sxs_profile_read(exp_path);
	
	double mean_radius = (A->rm * pa1->natoms + B->rm * pa2->natoms) / 
	                     (pa1->natoms + pa2->natoms);
	                     
	struct sxs_opt_params* params = sxs_opt_params_create(exp_profile, qvals, qnum, mean_radius);
	
	mol_prms_free(prms);
	free(exp_profile->qvals);
	sxs_profile_free(exp_profile);
	mol_atom_group_free(pa1);
	mol_atom_group_free(pa2);
	
//=================================================================
//=========================== correlate ===========================

	double* zvals = calloc(znum_total, sizeof(double));
	int  znum = 0;
	for (double zval = z_beg + myrank * z_step; 
	            zval < (z_end + 0.001); 
	            zval += nprocs * z_step) {
	            
		zvals[znum++] = zval;
	}
	
	if (myrank == 0) {
		SXS_PRINTF("\nCORRELATION STARTED\n\n");
	}
	
	// Compute scores
	sxs_fft_scores_full(A, B, params, 
	                    qvals, qnum, 
	                    zvals, znum, L);
	
	if (myrank == 0) {
		SXS_PRINTF("\nCORRELATION FINISHED\n\n");
	}
	
	sxs_opt_params_free(params);
	sxs_myfree(zvals);

#ifdef _MPI_
	if ((MPI_Finalize()) != MPI_SUCCESS)
	{
		fprintf (stderr, "MPI_Finalize failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
#endif

	return EXIT_SUCCESS;
}
