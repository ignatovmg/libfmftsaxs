#include "common.h"

#include "mol2/atom_group.h"
#include "mol2/pdb.h"
#include "mol2/prms.h"

#include "saxs_utils.h"
#include "form_factor_table.h"
#include "pdb2spf.h"
#include "profile.h"
#include "min_saxs.h"

void print_usage(char *exe_name)
{
	fprintf(stderr,
		  "Usage: %s MAPPING_PRM ATOMPRM REC_PATH LIG_PATH "
		  "EXP_PATH L PROFILE_PATH\n",
		  basename(exe_name));

	exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
	char* map_path = argv[1];
	char* prm_path = argv[2];
	char* rec_path = argv[3];
	char* lig_path = argv[4];
	char* exp_path = argv[5];
	int   L        = atoi(argv[6]);
	char* out_path = argv[7];
	
	if (argc != 8) { print_usage("single_saxs_fit"); }
	
	int   qnum = 50;
	double* qvals = sxs_mkarray(0.0, 0.5, qnum); 

	SXS_PRINTF("Reading parameters ...\n");
	struct mol_prms *prms = mol_prms_read(prm_path);
	
	SXS_PRINTF("Reading form-factors ...\n");
	struct saxs_form_factor_table *ff_table = default_ff_table(map_path);

	SXS_PRINTF("Reading receptor ...\n");
	struct mol_atom_group *rec = mol_read_pdb(rec_path);
	mol_atom_group_add_prms(rec, prms);
	
	SXS_PRINTF("Reading ligand ...\n");
	struct mol_atom_group *lig = mol_read_pdb(lig_path);
	mol_atom_group_add_prms(lig, prms);
	
	SXS_PRINTF("Reading experiment ...\n");
	struct sxs_profile* exp_profile = sxs_profile_read(exp_path);
	
	struct mol_atom_group *join = mol_atom_group_join(rec, lig);
	double join_rad = mol_atom_group_average_radius(join);
	struct sxs_opt_params* params = sxs_opt_params_create(exp_profile, qvals, qnum, join_rad);

	struct mol_vector3 com;
	centroid(&com, join);
	MOL_VEC_MULT_SCALAR(com, com, -1.0);
	mol_atom_group_translate(join, &com);

	SXS_PRINTF("Computing coefficients ...\n");
	struct sxs_profile* profile = sxs_profile_create(qvals, qnum, 1);
	struct sxs_spf_full *spf = atom_grp2spf(join, ff_table, qvals, qnum, L, 1);
	
	SXS_PRINTF("Fitting parameters ...\n\n");
	sxs_spf2fitted_profile(profile, spf, params);
	
	printf("Score: %.3f\nc1   : %.3f\nc2   : %.3f\n\n", 
	       profile->score, profile->c1, profile->c2);
	       
	sxs_profile_write(out_path, profile);
	
	SXS_PRINTF("Profile is written into %s\n", out_path);
	
	mol_atom_group_free(join);
	mol_atom_group_free(rec);
	mol_atom_group_free(lig);
	sxs_profile_free(profile);
	sxs_spf_full_free(spf);
	
	return EXIT_SUCCESS;
}
