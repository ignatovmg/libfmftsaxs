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
          "Usage: %s MAPPING_PRM ATOMPRM FTFILE RMFILE REC LIG "
          "SAXS_PROFILE OUTPUT L\n",
          basename(exe_name));
  
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	char* map_path = argv[1];
	char* prm_path = argv[2];
	char* ft__path = argv[3];
	char* rm__path = argv[4];
	char* rec_path = argv[5];
	char* lig_path = argv[6];
	char* exp_path = argv[7];
	int   L = atoi(argv[8]);
	char* out_path = argv[9];
	
	double z_beg  = 1.0;
	double z_end  = 80.0;
	double z_step = 1.0;
	
	if (argc != 10) { print_usage("correlate"); }
	
	
	int myrank = 0;
	int nprocs = 1;
	char eul_path[] = "eulertmp";
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
	double* qvals = mkarray(0.0, QMAX, qnum); 

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
	
	struct mol_vector3 ref_lig;
	MOL_VEC_SUB(ref_lig, com, coe);
	MOL_VEC_MULT_SCALAR(ref_lig, ref_lig, -1.0);
	sxs_ft_file2euler_file (eul_path, ft__path, rm__path, &ref_lig);
	
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
//=========================== massha ==============================  	  
	int* rec_num_list = calloc(nprocs, sizeof(int));
	
	double* zvals = calloc(znum_total, sizeof(double));
	int  znum = 0;
	for (double zval = z_beg + myrank * z_step; 
	            zval < (z_end + 0.001); 
	            zval += nprocs * z_step) {
	            
		zvals[znum++] = zval;
	}

	int*       ft_list_glb;
	double*   sx_score_glb;
	double*    c1_list_glb;
	double*    c2_list_glb;
	int ft_num_loc = 0;    
	int ft_num_glb = 0; 
	int ft_id;  
	int*    ft_list_loc;   
	double* sx_score_loc;
	double* c1_list_loc;
	double* c2_list_loc;
	
	int* sx_index_loc; 
	int* order_loc;
	int* order_glb;
	
	if (myrank == 0) { 
		SXS_PRINTF("Converting FT and RM files into Euler coordinates ...\n");
	}
	
	FILE* euler_file = fopen(eul_path, "r");
	
	double line[5], z;
	while(fscanf(euler_file, "%d %lf %lf %lf %lf %lf %lf", &ft_id, &z, 
	             &line[0], &line[1], &line[2], &line[3], &line[4]) != EOF)
	{
		for (int j = 0; j < znum; j++) {
			if (zvals[j] > z - 0.001 && zvals[j] < z + 0.001) {
				rec_num_list[myrank]++;
			}
		}
	}
	rewind(euler_file);
	
#ifdef _MPI_
	for (int i = 0; i < nprocs; i++) {
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&rec_num_list[i], 1, MPI_INT, i, MPI_COMM_WORLD);
	}
#endif
	
	ft_num_loc = rec_num_list[myrank];
	for (int i = 0; i < nprocs; i++) { 
		ft_num_glb += rec_num_list[i]; 
	}
	
	sx_index_loc = calloc(ft_num_loc, sizeof(int));
	sx_score_loc = calloc(ft_num_loc, sizeof(double));
	ft_list_loc = calloc(ft_num_loc, sizeof(int));
	c1_list_loc = calloc(ft_num_loc, sizeof(double));
	c2_list_loc = calloc(ft_num_loc, sizeof(double));
	order_loc = calloc(ft_num_loc, sizeof(int));

	double b_step = M_PI / L;
	double a_step = 2.0 * M_PI / (2*L+1);
	int counter = 0, iline = 0;
	int tmp, tmp_id;


	// Read the file with euler coordinates of the conformations
	// and find the index of the closest grid point for each of
	// them
	
	if (myrank == 0) { SXS_PRINTF("Reading Euler coordinates ...\n"); }
	
	struct sxs_euler euler;
	while(fscanf(euler_file, "%d %lf %lf %lf %lf %lf %lf", &ft_id, 
	             &euler.z,  &euler.b1, &euler.g1, &euler.a2, 
	             &euler.b2, &euler.g2) != EOF) {

		for (int j = 0; j < znum; j++) {
			if (zvals[j] > euler.z - 0.001 && zvals[j] < euler.z + 0.001) {
			
				euler.a2 = 2*M_PI - euler.a2;
				euler.g2 = 2*M_PI - euler.g2;
				//euler.g1 = 2*M_PI - euler.g2;
	
				tmp_id = j * nbeta;

				tmp = (int)(round(euler.b1 / b_step));
				tmp_id = (tmp_id + tmp) * nbeta;

				tmp = (int)(round(euler.b2 / b_step));
				tmp_id = (tmp_id + tmp) * (2*L+1);

				tmp = (int)(round(euler.a2 / a_step));
				tmp_id = (tmp_id + tmp) * (2*L+1);

				tmp = (int)(round(euler.g1 / a_step));
				tmp_id = (tmp_id + tmp) * (2*L+1);

				tmp = (int)(round(euler.g2 / a_step));
				tmp_id = tmp_id + tmp;

				sx_index_loc[counter] = tmp_id;
				ft_list_loc[counter] = ft_id;
				order_loc[counter] = iline;
				
				counter++;
			}
		}
		
		iline++;
	}

	fclose(euler_file);

	if (myrank == 0)
	{
		// Here the scores computed by the other processes will
		// be stored
	
		sx_score_glb  = calloc(ft_num_glb, sizeof(double));
		ft_list_glb   = calloc(ft_num_glb, sizeof(int));
		c1_list_glb   = calloc(ft_num_glb, sizeof(double));
		c2_list_glb   = calloc(ft_num_glb, sizeof(double));
		order_glb     = calloc(ft_num_glb, sizeof(int));
	}
	
	if (myrank == 0) {
		SXS_PRINTF("\nCORRELATION STARTED\n\n");
	}
	
	// Compute scores
	compute_saxs_scores(sx_score_loc, 
	                    c1_list_loc, 
	                    c2_list_loc,
	                    sx_index_loc, 
	                    ft_num_loc, 
	                    A, B, 
	                    params, 
	                    qvals, qnum, 
	                    zvals, znum, 
	                    L, 1);
		
#ifndef _MPI_
	endtime = (double)clock() / CLOCKS_PER_SEC;
	
	sx_score_glb = sx_score_loc;
	ft_list_glb = ft_list_loc;
	c1_list_glb = c1_list_loc;
	c2_list_glb = c2_list_loc;
	order_glb   = order_loc;
#else
	// Send computed values to the master
	
	// Fill arrays for MPI_Gatherv()
	int displacements_array[nprocs];
	int proc_out_size_array[nprocs];
	
	memset(displacements_array, 0, nprocs*sizeof(int));
	for (int i = 0; i < nprocs; i++) {
		if (i > 0) {
			displacements_array[i] = displacements_array[i - 1] + rec_num_list[i - 1];
		}
	}
	
	// collect ft indices
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(ft_list_loc, ft_num_loc, MPI_INT, ft_list_glb, rec_num_list, displacements_array, MPI_INT, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
	
	sxs_myfree(ft_list_loc);
	
	// collect order
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(order_loc, ft_num_loc, MPI_INT, order_glb, rec_num_list, displacements_array, MPI_INT, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
	
	sxs_myfree(order_loc);
	
	// collect scores
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(sx_score_loc, ft_num_loc, MPI_DOUBLE, sx_score_glb, rec_num_list, displacements_array, MPI_DOUBLE, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
	
	sxs_myfree(sx_score_loc);

    // collect c1
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(c1_list_loc, ft_num_loc, MPI_DOUBLE, c1_list_glb, rec_num_list, displacements_array, MPI_DOUBLE, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}

	sxs_myfree(c1_list_loc);

    // collect c2
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(c2_list_loc, ft_num_loc, MPI_DOUBLE, c2_list_glb, rec_num_list, displacements_array, MPI_DOUBLE, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
	
	sxs_myfree(c2_list_loc);

	MPI_Barrier(MPI_COMM_WORLD);

	endtime = MPI_Wtime();
#endif


	// Master records the scores into out_path
	if (myrank == 0) {
	
		printf("\nTime passed: %.3f\n", endtime - starttime);
		printf("Writing results to %s\n", out_path);
		FILE* stream;
		
		stream = fopen(out_path, "w");
		for (int i = 0; i < ft_num_glb; i++) {
			fprintf(stream, "%6.d\t%d\t%.3lf\t%.3lf\t%.3lf\n", order_glb[i], 
			                                              ft_list_glb[i], 
			                                             sx_score_glb[i], 
			                                             c1_list_glb[i], 
			                                             c2_list_glb[i]);
		}

		fclose(stream);
		
		sxs_myfree(sx_score_glb);
		sxs_myfree(ft_list_glb);
		sxs_myfree(c1_list_glb);
		sxs_myfree(c2_list_glb);
		sxs_myfree(order_glb);
		
		printf("\nCorrelation finished\n");
	}
	
	sxs_opt_params_free(params);
	sxs_myfree(rec_num_list);
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
