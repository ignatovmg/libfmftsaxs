
#include "common.h"
#include <time.h>

#include "massha2.h"
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
          "SAXS_PROFILE OUTPUT L ZSTART ZEND ZSTEP\n",
          basename(exe_name));
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
	char* out_path = argv[8];
	
	int L = atoi(argv[9]);
	double z_beg  = atof(argv[10]);
	double z_end  = atof(argv[11]);
	double z_step = atof(argv[12]);
	
	if (argc != 13) {
		print_usage("fitmassha2_new");
	}
	
	
	int myrank = 0;
	int nprocs = 1;
	char eul_path[] = "eutmp";
	
	int nbeta = L + 1;
	int znum_total = (int)round((z_end - z_beg) / z_step) + 1;
	double overalltime = 0.0;
	
	int qnum = QNUM;                                      // number of points, in which saxs curve is evaluated
	double* qvals = mkarray(0.0, QMAX, qnum); 

	struct mol_prms *prms = mol_prms_read(prm_path);
	struct saxs_form_factor_table *ff_table = default_ff_table(map_path);
	
//=================================================================
//=========================== rec =================================
	struct mol_atom_group *pa1 = mol_read_pdb(rec_path);
	mol_atom_group_add_prms(pa1, prms);
	
	struct mol_vector3 coe;
	center_of_extrema(&coe, pa1);
	MOL_VEC_MULT_SCALAR(coe, coe, -1.0);
	mol_atom_group_translate(pa1, &coe);

	struct sxs_spf_full* A = atom_grp2spf(pa1, ff_table, qvals, qnum, L, 1);

	printf("\n_______________________\nFirst protein finished\n_______________________\n\n");  

//=================================================================
//=========================== lig =================================
	struct mol_atom_group *pa2 = mol_read_pdb(lig_path);
	mol_atom_group_add_prms(pa2, prms);
	
	struct mol_vector3 com;
	centroid(&com, pa2);
	MOL_VEC_MULT_SCALAR(com, com, -1.0);
	mol_atom_group_translate(pa2, &com);
	
	struct sxs_spf_full* B = atom_grp2spf(pa2, ff_table, qvals, qnum, L, 1);
	
	printf("\n_______________________\nSecond protein finished\n_______________________\n\n");   
	
	struct mol_vector3 ref_lig;
	MOL_VEC_SUB(ref_lig, com, coe);
	MOL_VEC_MULT_SCALAR(ref_lig, ref_lig, -1.0);
	saxs_ft_file2euler_file (eul_path, ft__path, rm__path, &ref_lig);
	
//=================================================================
//=================== experiment ==================================
	struct sxs_profile* exp_profile = sxs_profile_read(exp_path);
	
	double mean_radius = (A->rm * pa1->natoms + B->rm * pa2->natoms) / 
	                     (pa1->natoms + pa2->natoms);
	                     
	struct sxs_opt_params* params = sxs_opt_params_create(exp_profile, qvals, qnum, mean_radius);
	
	mol_prms_free(prms);
	sxs_profile_free(exp_profile);

	printf("\n_______________________\n'Experimental' protein finished\n_______________________\n\n"); 

//=================================================================
//=========================== massha ==============================  	  


	printf("\nCORRELATION STARTED\n\n");
	clock_t time = clock();
	
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
#endif
	
	int znum_proc;
	int znum_list[nprocs];
	
	for (int i = 0; i < nprocs; i++) {
		znum_list[i] = znum_total / nprocs;
		
		if (i < znum_total % nprocs) {
			znum_list[i]++;
		}
	}

	double* z_beg_list = calloc(nprocs, sizeof(double));
	double* z_end_list = calloc(nprocs, sizeof(double));
	
	z_beg_list[0] = z_beg;
	z_end_list[0] = z_beg + (znum_list[0] - 1) * z_step;
	for (int i = 1; i < nprocs; i++) {
		z_beg_list[i] = z_end_list[i-1] + z_step;
		z_end_list[i] = z_beg_list[i] + (znum_list[i] - 1) * z_step;
	}
	
	znum_proc = znum_list[myrank];
	z_beg = z_beg_list[myrank];
	z_end = z_end_list[myrank];
	

	int*       ft_list_glb;
	double*   sx_score_glb;
	double*    c1_list_glb;
	double*    c2_list_glb;
	int ft_num_loc = 0;    
	int ft_num_glb = 0; 
	int ft_id;
	int* rec_path_num_list;    
	int*    ft_list_loc;   
	double* sx_score_loc;
	double* c1_list_loc;
	double* c2_list_loc;
	
	int* saxs_index_loc; 
	
	FILE* euler_file = fopen(eul_path, "r");
	rec_path_num_list = calloc(nprocs, sizeof(int));
	
	double line[5], z;
	while(fscanf(euler_file, "%d %lf %lf %lf %lf %lf %lf", &ft_id, &z, 
	             &line[0], &line[1], &line[2], &line[3], &line[4]) != EOF)
	{
		for (int i = 0; i < nprocs; i++) {
			if (z >= z_beg_list[i] && z <= z_end_list[i]) {
				rec_path_num_list[i]++;
				ft_num_glb++;
			}
		}
	}
	rewind(euler_file);
	ft_num_loc = rec_path_num_list[myrank];
	
	saxs_index_loc = calloc(ft_num_loc, sizeof(int));
	sx_score_loc = calloc(ft_num_loc, sizeof(double));
	ft_list_loc = calloc(ft_num_loc, sizeof(int));
	c1_list_loc = calloc(ft_num_loc, sizeof(double));
	c2_list_loc = calloc(ft_num_loc, sizeof(double));
	
	double b_step = M_PI / L;
	double a_step = 2.0 * M_PI / (2*L+1);
	int counter = 0;
	int tmp, tmp_id;
	
	// convert euler coordinates into complex ids in saxs array
	struct saxs_euler euler;
	while(fscanf(euler_file, "%d %lf %lf %lf %lf %lf %lf", &ft_id, 
	             &euler.z,  &euler.b1, &euler.g1, &euler.a2, 
	             &euler.b2, &euler.g2) != EOF) {

		if (euler.z < z_beg || euler.z > z_end) {
			continue;
		}
		
		euler.a2 = 2*M_PI - euler.a2;
		euler.g2 = 2*M_PI - euler.g2;
		//euler.g1 = 2*M_PI - euler.g2;
	
        tmp_id = (int)(round((euler.z - z_beg) / z_step)) * nbeta;

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

        saxs_index_loc[counter] = tmp_id;
        ft_list_loc[counter++] = ft_id;
	}
	
	fclose(euler_file);
	
	if (myrank == 0)
	{
		sx_score_glb  = calloc(ft_num_glb, sizeof(double));
		ft_list_glb = calloc(ft_num_glb, sizeof(int));
		c1_list_glb = calloc(ft_num_glb, sizeof(double));
		c2_list_glb = calloc(ft_num_glb, sizeof(double));
	}
	
	compute_saxs_scores(
		A, B, qvals, qnum, 
		params,
		z_beg, z_step, 
		znum_proc,
		sx_score_loc,
		saxs_index_loc,
		c1_list_loc,
		c2_list_loc,
		ft_num_loc, L);
		
#ifndef _MPI_
	time = clock() - time;
	overalltime += (double)time / CLOCKS_PER_SEC;
	
	sx_score_glb = sx_score_loc;
	ft_list_glb = ft_list_loc;
	c1_list_glb = c1_list_loc;
	c2_list_glb = c2_list_loc;
#else
	//printf("For process %i calculation finished -> time passed: %f\n", myrank, processtime);
	
	// Fill arrays for MPI_Gatherv()
	int displacements_array[nprocs];
	int proc_out_size_array[nprocs];
	
	memset(displacements_array, 0, nprocs*sizeof(int));
	for (int i = 0; i < nprocs; i++) {
		if (i > 0) {
			displacements_array[i] = displacements_array[i - 1] + rec_path_num_list[i - 1];
		}
	}
	
	// collect ft indices
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(ft_list_loc, ft_num_loc, MPI_INT, ft_list_glb, rec_path_num_list, displacements_array, MPI_INT, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
	
	// collect scores
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(sx_score_loc, ft_num_loc, MPI_DOUBLE, sx_score_glb, rec_path_num_list, displacements_array, MPI_DOUBLE, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}

    // collect c1
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(c1_list_loc, ft_num_loc, MPI_DOUBLE, c1_list_glb, rec_path_num_list, displacements_array, MPI_DOUBLE, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}

    // collect c2
	MPI_Barrier(MPI_COMM_WORLD);
	
	if ((MPI_Gatherv(c2_list_loc, ft_num_loc, MPI_DOUBLE, c2_list_glb, rec_path_num_list, displacements_array, MPI_DOUBLE, 0, MPI_COMM_WORLD)) != MPI_SUCCESS) {
		fprintf (stderr, "MPI_Gatherv failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	free(rec_path_num_list);
	free(ft_list_loc);
	free(sx_score_loc);
	free(saxs_index_loc);
	
	overalltime = MPI_Wtime() - overalltime;
	
#endif

	if (myrank == 0) {
	
		printf("Overall time passed: %f\n", overalltime);

		printf("Writing results...\n");
		FILE* stream;
		
		stream = fopen(out_path, "w");
		for (int i = 0; i < ft_num_glb; i++) {
			fprintf(stream, "%d\t%.3lf\t%.3lf\t%.3lf\n", ft_list_glb[i], 
			                                             sx_score_glb[i], 
			                                             c1_list_glb[i], 
			                                             c2_list_glb[i]);
		}

		fclose(stream);
		
		free(sx_score_glb);
		free(ft_list_glb);
		free(c1_list_glb);
		free(c2_list_glb);
		
		printf("\n_______________________\nSAXS finished\n_______________________\n\n");
	}
	
	sxs_opt_params_free(params);

#ifdef _MPI_
	if ((MPI_Finalize()) != MPI_SUCCESS)
	{
		fprintf (stderr, "MPI_Finalize failed\n");
		MPI_Abort (MPI_COMM_WORLD, 1);
	}
#endif

	return 0;
}
