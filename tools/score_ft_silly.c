#include "common.h"

#include "fftsaxs.h"
#include "index.h"
#include "saxs_utils.h"

#include "mol2/atom_group.h"
#include "mol2/pdb.h"
#include "mol2/prms.h"


void rotate(struct mol_atom_group *pa, double x0, double y0, double z0, double fi)
{
    int natoms = pa->natoms;
    double r = sqrt(x0*x0+y0*y0+z0*z0);
    double x = x0/r;
    double y = y0/r;
    double z = z0/r;
    double r_matrix[3][3];
    r_matrix[0][0] = (cos(fi) + (1 - cos(fi)) * x * x); 
    r_matrix[0][1] = ((1 - cos(fi)) * x * y - sin(fi) * z); 
    r_matrix[0][2] = ((1 - cos(fi)) * x * z + sin(fi) * y);
    r_matrix[1][0] = ((1 - cos(fi)) * x * y + sin(fi) * z); 
    r_matrix[1][1] = (cos(fi) + (1 - cos(fi)) * y * y); 
    r_matrix[1][2] = ((1 - cos(fi)) * z * y - sin(fi) * x);
    r_matrix[2][0] = ((1 - cos(fi)) * x * z - sin(fi) * y); 
    r_matrix[2][1] = ((1 - cos(fi)) * z * y + sin(fi) * x); 
    r_matrix[2][2] = (cos(fi) + (1 - cos(fi)) * z * z);
    
    for (int i = 0; i < natoms; i++)
    {       
        double new_x = r_matrix[0][0] * pa->coords[i].X + 
        r_matrix[0][1] * pa->coords[i].Y + 
        r_matrix[0][2] * pa->coords[i].Z;
        
        double new_y = r_matrix[1][0] * pa->coords[i].X + 
        r_matrix[1][1] * pa->coords[i].Y + 
        r_matrix[1][2] * pa->coords[i].Z;
        
        double new_z = r_matrix[2][0] * pa->coords[i].X + 
        r_matrix[2][1] * pa->coords[i].Y + 
        r_matrix[2][2] * pa->coords[i].Z;
        
        pa->coords[i].X = new_x;
        pa->coords[i].Y = new_y;
        pa->coords[i].Z = new_z;
    }
}

void euler_rotate_zyz(struct mol_atom_group *pa, double alpha, double beta, double gamma)
{
	rotate(pa, 0.0, 0.0, 1.0, gamma);
	rotate(pa, 0.0, 1.0, 0.0, beta);
	rotate(pa, 0.0, 0.0, 1.0, alpha);
}

void print_usage(char *exe_name)
{
	fprintf(stderr,
		  "Usage: %s MAPPING_PRM ATOMPRM FT_PATH RM_PATH "
		  " REC_PATH LIG_PATH REF_PROFILE L\n",
		  basename(exe_name));

	exit(EXIT_FAILURE);
}


int main(int argc, char** argv)
{
	char* map_path = argv[1];
	char* prm_path = argv[2];
	char* ft_path  = argv[3];
	char* rm_path  = argv[4];
	char* rec_path = argv[5];
	char* lig_path = argv[6];
	char* exp_path = argv[7];
	int   L = atoi(argv[8]);
	char* out_path = argv[9];
	int   qnum = 50;
	
	if (argc != 10) { print_usage("score_ft_silly"); }
	
	char eul_path[] = "euler_list";
	double* qvals = mkarray(0.0, 0.5, qnum); 

	struct mol_prms *prms = mol_prms_read(prm_path);
	struct saxs_form_factor_table *ff_table = default_ff_table(map_path);

	// read rec and move to coe
	struct mol_atom_group *rec = mol_read_pdb(rec_path);
	mol_atom_group_add_prms(rec, prms);
	struct mol_atom_group *rec_moved = mol_read_pdb(rec_path);
	mol_atom_group_add_prms(rec_moved, prms);
	
	struct mol_vector3 coe;
	center_of_extrema(&coe, rec);
	MOL_VEC_MULT_SCALAR(coe, coe, -1.0);
	mol_atom_group_translate(rec, &coe);
	mol_atom_group_translate(rec_moved, &coe);
	
	// read lig and move to com
	struct mol_atom_group *lig = mol_read_pdb(lig_path);
	mol_atom_group_add_prms(lig, prms);
	struct mol_atom_group *lig_moved = mol_read_pdb(lig_path);
	mol_atom_group_add_prms(lig_moved, prms);
	
	struct mol_vector3 com;
	centroid(&com, lig);
	MOL_VEC_MULT_SCALAR(com, com, -1.0);
	mol_atom_group_translate(lig, &com);
	mol_atom_group_translate(lig_moved, &com);
	
	// reference ligand position
	struct mol_vector3 ref_lig;
	MOL_VEC_SUB(ref_lig, com, coe);
	MOL_VEC_MULT_SCALAR(ref_lig, ref_lig, -1.0);
	
	// read experiment
	struct sxs_profile* exp_profile = sxs_profile_read(exp_path);
	
	// set minimizer params
	double rec_rad = mol_atom_group_average_radius(rec);
	double lig_rad = mol_atom_group_average_radius(lig);
	double join_rad = (rec_rad * rec->natoms + lig_rad * lig->natoms) / (rec->natoms + lig->natoms);
	struct sxs_opt_params* params = sxs_opt_params_create(exp_profile, qvals, qnum, join_rad);
	
	struct sxs_profile* profile = sxs_profile_create(qvals, qnum, 1);
	
	// record and open euler file
	sxs_ft_file2euler_file(eul_path, ft_path, rm_path, &ref_lig);
	FILE* euler = sxs_myfopen(eul_path, "r");
	
	int nread, eu_id;
	double z, b1, g1, a2, b2, g2;
	struct mol_matrix3 rec_rm;
	struct mol_matrix3 lig_rm;
	struct mol_vector3 rec_tv = {0.0, 0.0, 0.0};
	struct mol_vector3 lig_tv = {0.0, 0.0, 0.0};
	
	struct sxs_spf_full *coefs = sxs_spf_full_create(L, qnum);
	
	// join groups and compute solvent accessibility for each atom
	struct mol_atom_group *join = mol_atom_group_join(rec_moved, lig_moved);
	double* saxs_sa = malloc(join->natoms * sizeof(double));
	faccs(saxs_sa, join, 1.4);
	mol_atom_group_free(join);
	
	FILE* out = fopen(out_path, "w");
	
	clock_t time1, time2;
	time2 = clock();
	int iline = 0;
	
	printf("Output will be written to %s\n", out_path);
	
	while (fscanf(euler, "%d %lf %lf %lf %lf %lf %lf", &eu_id, &z, 
	             &b1, &g1, &a2, &b2, &g2) != EOF) {
	        
	    printf("Processing line %i\n", iline++);
	        
	    SXS_PRINTF("euid = %i\n", eu_id);
		time1 = clock();
	   
	    // arrange and join atom groups
		lig_tv.Z = z;          
		sxs_fill_active_rotation_matrix(&rec_rm, 0.0, b1, g1);
		sxs_fill_active_rotation_matrix(&lig_rm,  a2, b2, g2);
		mol_atom_group_move_in_copy(lig, lig_moved, &lig_rm, &lig_tv);
		mol_atom_group_move_in_copy(rec, rec_moved, &rec_rm, &rec_tv);
		
		SXS_PRINTF("\tmove: %.4f\n", (double)(clock()-time1) / CLOCKS_PER_SEC);     
		time1 = clock();
		
		join = mol_atom_group_join(rec_moved, lig_moved);
		
		SXS_PRINTF("\tjoin: %.4f\n", (double)(clock()-time1) / CLOCKS_PER_SEC);     
		time1 = clock();
		
		// move to the com
		struct mol_vector3 joincom;
		centroid(&joincom, join);
		MOL_VEC_MULT_SCALAR(joincom, joincom, -1.0);
		mol_atom_group_translate(join, &joincom);
		
		SXS_PRINTF("\tcntr: %.4f\n", (double)(clock()-time1) / CLOCKS_PER_SEC);     
		time1 = clock();
		
		// compute SPF decomposition (see pdb2spf.h)
		atom_grp2spf_inplace(coefs, join, ff_table, qvals, qnum, L, saxs_sa);
		
		SXS_PRINTF("\tampl: %.4f\n", (double)(clock()-time1) / CLOCKS_PER_SEC);     
		time1 = clock();
		
		//sxs_profile_from_spf(profile, coefs, 1.0, 0.0);
		
		// assemble and minimize profile
		spf2fitted_profile(profile, coefs, params);
		
		SXS_PRINTF("\tfit : %.4f\n", (double)(clock()-time1) / CLOCKS_PER_SEC);     
		time1 = clock();
		
		SXS_PRINTF("\tchi :%d\t%.3f\t%.3f\t%.3f\n", eu_id, profile->score, profile->c1, profile->c2);
		fprintf(out, "%d\t%.3f\t%.3f\t%.3f\n", eu_id, profile->score, profile->c1, profile->c2);

		mol_atom_group_free(join);
		
		SXS_PRINTF("\tfree: %.4f\n\n", (double)(clock()-time1) / CLOCKS_PER_SEC);     
		time1 = clock();

	}
	
	printf("Total time elapsed: %.4f\n", (double)(clock()-time2) / CLOCKS_PER_SEC);     

	fclose(out);
	fclose(euler);
	
	sxs_myfree(saxs_sa);
	mol_prms_free(prms);
	sxs_spf_full_free(coefs);
	sxs_profile_free(exp_profile);
	sxs_profile_free(profile);
	sxs_opt_params_free(params);
	mol_atom_group_free(rec);
	mol_atom_group_free(lig);
	free(ff_table);	

	return 0;
}
