#include "common.h"

#include "mol2/atom_group.h"
#include "mol2/pdb.h"
#include "mol2/prms.h"

#include "saxs_utils.h"
#include "form_factor_table.h"
#include "pdb2spf.h"
#include "mpapi2.h"
#include "min_saxs.h"

#include <gmp.h>

int main(int argc, char** argv)
{
	/*char pdb[] = "data/4g9s_r_native.pdb";
	struct mol_atom_group* ag = mol_read_pdb(pdb);
	struct mol_prms *prms = mol_prms_read("mol-prms/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec");
	mol_atom_group_add_prms(ag, prms);
	
	struct mol_vector3 center;
	center_of_extrema(&center, ag);
	MOL_VEC_MULT_SCALAR(center, center, -1.0);
	mol_atom_group_translate(ag, &center);
	
	struct saxs_form_factor_table *ff_table = default_ff_table("mol-prms/pdb_formfactor_mapping_clean.prm");
	
	printf("%f\n", ff_table->factors[s_OH2].zero_ff);
	
	int L = 15;
	int qnum = 50;
	double* qvals = mkarray(0.0, 0.5, qnum);
	
	mpf_t fp, rop;
	mpf_init(rop);
	mpf_init_set_d(fp, 345.0);
	mpf_out_str(stdout, 10, 30, fp);
	printf("\n");
	mpf_pow_ui(rop, fp, 1000);
	mpf_out_str(stdout, 10, 3, rop);
	printf("\n");
	
	mpf_t* fpar = calloc(10, sizeof(mpf_t));
	for (int i = 0; i < 10; i++) {
		mpf_init(fpar[i]);
		mpf_set_d(fpar[i], 10.0);
		mpf_pow_ui(fpar[i], fpar[i], i);
		mpf_out_str(stdout, 10, 3, fpar[i]);
		printf("\n");
	}
	
	struct sxs_spf_full* spf = atom_grp2spf(ag, ff_table, qvals, qnum, L, 1);
	sxs_spf_full_write("my_spf.txt", spf);
	
	struct sxs_profile* profile = sxs_profile_create(qvals, qnum, 0);
	sxs_profile_from_spf(profile, spf, 1.0, 1.0);
	sxs_profile_write("my_profile.txt", profile);*/
	
	char mapping_file_path[] = "mol-prms/pdb_formfactor_mapping_clean.prm";          // for finding atoms in param file
	char prms_file_path[]    = "mol-prms/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec";          // param file
	char dimer[]          = "data/4g9s_native_dimer.pdb";
	char exp_path[]		= "data/ref_saxs";
	int  L              = 15;
	int  qnum = 50;
	double* qvals = mkarray(0.0, 0.5, qnum); 

	struct mol_prms *prms = mol_prms_read(prms_file_path);
	struct saxs_form_factor_table *ff_table = default_ff_table(mapping_file_path);

	struct mol_atom_group *pa = mol_read_pdb(dimer);
	mol_atom_group_add_prms(pa, prms);

	struct mol_vector3 coe;
	center_of_extrema(&coe, pa);
	MOL_VEC_MULT_SCALAR(coe, coe, -1.0);
	mol_atom_group_translate(pa, &coe);

	struct sxs_spf_full* A = atom_grp2spf(pa, ff_table, qvals, qnum, L, 1);
	
	struct sxs_profile* exp_profile = sxs_profile_read(exp_path);
	struct sxs_opt_params* params = sxs_opt_params_create(exp_profile, qvals, qnum, A->rm);

	struct sxs_profile* profile = sxs_profile_create(qvals, qnum, 1);
	sxs_profile_from_spf(profile, A, 1.0, 0.0);
	sxs_profile_write("unfitted_profile", profile);
	
	spf2fitted_profile(profile, A, params);
	sxs_profile_write("fitted_profile", profile);

	printf("%f\n", profile->score);

	mol_prms_free(prms);
	sxs_profile_free(exp_profile);
	sxs_profile_free(profile);
	sxs_opt_params_free(params);
	mol_atom_group_free(pa);
	free(ff_table);
	
	return 0;
}
