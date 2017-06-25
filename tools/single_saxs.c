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
		  "C1 C2 L PROFILE_PATH\n",
		  basename(exe_name));

	exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
	char* map_path = argv[1];
	char* prm_path = argv[2];
	char* rec_path = argv[3];
	char* lig_path = argv[4];
	double c1 = atof(argv[5]);
	double c2 = atof(argv[6]);
	int    L  = atoi(argv[7]);
	char* out_path = argv[8];
	
	if (argc != 9) { print_usage("single_saxs"); }
	
	int   qnum = 50;
	double* qvals = mkarray(0.0, 0.5, qnum); 

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
	
	struct mol_atom_group *join = mol_atom_group_join(rec, lig);

	struct mol_vector3 com;
	centroid(&com, join);
	MOL_VEC_MULT_SCALAR(com, com, -1.0);
	mol_atom_group_translate(join, &com);

	SXS_PRINTF("Computing coefficients ...\n");
	struct sxs_profile* profile = sxs_profile_create(qvals, qnum, 1);
	struct sxs_spf_full *spf = atom_grp2spf(join, ff_table, qvals, qnum, L, 1);
	sxs_profile_from_spf(profile, spf, c1, c2);

	sxs_profile_write(out_path, profile);
	
	SXS_PRINTF("Profile is written into %s\n", out_path);
	
	mol_atom_group_free(join);
	mol_atom_group_free(rec);
	mol_atom_group_free(lig);
	sxs_profile_free(profile);
	sxs_spf_full_free(spf);
	
	return EXIT_SUCCESS;
}
