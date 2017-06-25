#include "common.h"

#include "mol2/atom_group.h"
#include "mol2/pdb.h"
#include "mol2/prms.h"
#include "pdb2spf.h"
#include "profile.h"

#include "fftsaxs.h"
#include "index.h"

#ifdef _MPI_
#include <mpi.h>
#endif

#include <check.h>


const double tol1 = 1e-3;
const double tol2 = 1e-6;
const double tol3 = 1e-9;

static inline void double_eq_rel_tol(double val1, double val2, double tol)
{
	if (fabs(val1) > 0.0 && fabs(val2) > 0.0) {
		ck_assert_double_le_tol(fabs((val2 - val1) / val2), tol, tol3);
	} else {
		ck_assert_double_eq_tol(val2, val1, tol2);
	}
}

static void cmp_profiles(struct sxs_profile* cur_profile, struct sxs_profile* ref_profile)
{
	ck_assert_int_eq(cur_profile->qnum, ref_profile->qnum);
	int qnum = cur_profile->qnum;

	double val1, val2;
	for (int q = 0; q < qnum; q++) {
		//printf("% .4f ^ % .4f, % .3f ^ % .3f, % .3f ^ % .3f\n", 
		//         cur_profile->qvals[q], ref_profile->qvals[q], 
		//         cur_profile->in[q],    ref_profile->in[q], 
		//         cur_profile->err[q],   ref_profile->err[q]);
		
		val1 = cur_profile->qvals[q];
		val2 = ref_profile->qvals[q];
		ck_assert_double_eq_tol(val1, val2, 10e-4);
		
		val1 = cur_profile->in[q];
		val2 = ref_profile->in[q];
		double_eq_rel_tol(val1, val2, tol2);
		
		val1 = cur_profile->err[q];
		val2 = ref_profile->err[q];
		double_eq_rel_tol(val1, val2, tol2);
	}
}

// Test cases
START_TEST(test_pdb2spf)
{
	struct mol_atom_group* ag = mol_read_pdb("4g9s_r_native.pdb");
	struct mol_prms *prms = mol_prms_read("atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec");
	
	ck_assert_ptr_nonnull(ag);
	ck_assert_ptr_nonnull(prms);
	
	mol_atom_group_add_prms(ag, prms);
	
	struct mol_vector3 center;
	center_of_extrema(&center, ag);
	MOL_VEC_MULT_SCALAR(center, center, -1.0);
	mol_atom_group_translate(ag, &center);

	struct saxs_form_factor_table *ff_table = default_ff_table("pdb_formfactor_mapping_clean.prm");
	ck_assert_ptr_nonnull(ff_table);
	
	int L = 15;
	int qnum = 50;
	double* qvals = mkarray(0.0, 0.5, qnum);
	ck_assert_ptr_nonnull(qvals);
	
	struct sxs_spf_full* spf_coef = atom_grp2spf(ag, ff_table, qvals, qnum, L, 1);
	ck_assert_ptr_nonnull(spf_coef);
	
	sxs_spf_full_write("../ref_spf", spf_coef);
	
	struct sxs_spf_full* spf_coef_ref = sxs_spf_full_read("ref_spf");
	ck_assert_ptr_nonnull(spf_coef_ref);
	
	double val1, val2;
	for (int q = 0; q < qnum; q++) {
		for (int i = 0; i < (L+1)*(L+1); i++) {
		
			val1 = spf_coef->V[q]->re[i];
			val2 = spf_coef_ref->V[q]->re[i];
			double_eq_rel_tol(val1, val2, tol1);
			
			val1 = spf_coef->V[q]->im[i];
			val2 = spf_coef_ref->V[q]->im[i];
			double_eq_rel_tol(val1, val2, tol1);
			
			val1 = spf_coef->D[q]->re[i];
			val2 = spf_coef_ref->D[q]->re[i];
			double_eq_rel_tol(val1, val2, tol1);
			
			val1 = spf_coef->D[q]->im[i];
			val2 = spf_coef_ref->D[q]->im[i];
			double_eq_rel_tol(val1, val2, tol1);
			
			val1 = spf_coef->W[q]->re[i];
			val2 = spf_coef_ref->W[q]->re[i];
			double_eq_rel_tol(val1, val2, tol1);
			
			val1 = spf_coef->W[q]->im[i];
			val2 = spf_coef_ref->W[q]->im[i];
			double_eq_rel_tol(val1, val2, tol1);
		}
	}
	
	//ff_table_destroy(ff_table);
	mol_atom_group_free(ag);
	mol_prms_free(prms);
	sxs_spf_full_free(spf_coef_ref);
	sxs_spf_full_free(spf_coef);
	sxs_myfree(qvals);
}
END_TEST

START_TEST(test_sxs_profile_from_spf)
{
	struct mol_atom_group* ag = mol_read_pdb("4g9s_r_native.pdb");
	struct mol_prms *prms = mol_prms_read("atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec");
	
	ck_assert_ptr_nonnull(ag);
	ck_assert_ptr_nonnull(prms);
	
	mol_atom_group_add_prms(ag, prms);
	
	struct mol_vector3 center;
	center_of_extrema(&center, ag);
	MOL_VEC_MULT_SCALAR(center, center, -1.0);
	mol_atom_group_translate(ag, &center);

	struct saxs_form_factor_table *ff_table = default_ff_table("pdb_formfactor_mapping_clean.prm");
	ck_assert_ptr_nonnull(ff_table);
	
	int L = 15;
	int qnum = 50;
	double* qvals = mkarray(0.0, 0.5, qnum);
	ck_assert_ptr_nonnull(qvals);
	
	struct sxs_spf_full* spf_coef = atom_grp2spf(ag, ff_table, qvals, qnum, L, 1);
	ck_assert_ptr_nonnull(spf_coef);

	struct sxs_profile* cur_profile = sxs_profile_create(qvals, qnum, 0);
	ck_assert_ptr_nonnull(cur_profile);
	
	sxs_profile_from_spf(cur_profile, spf_coef, 1.0, 1.0);
	
	sxs_profile_write("../ref_profile", cur_profile);
	
	struct sxs_profile* ref_profile = sxs_profile_read("ref_profile");
	ck_assert_ptr_nonnull(ref_profile);

	cmp_profiles(cur_profile, ref_profile);
	
	//ff_table_destroy(ff_table);
	sxs_profile_free(cur_profile);
	sxs_profile_free(ref_profile);
	sxs_spf_full_free(spf_coef);
	sxs_myfree(qvals);
	mol_atom_group_free(ag);
	mol_prms_free(prms);
}
END_TEST

START_TEST(score_conformations)
{
	char map_path[] = "pdb_formfactor_mapping_clean.prm"; 
	char prm_path[] = "atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec";
	char rec_path[] = "4g9s_r_native.pdb";
	char lig_path[] = "4g9s_l_moved.pdb"; 
	char exp_path[]	= "ref_saxs";
	char eul_path[] = "euler_coords.000.00";
	char ref_path[] = "ref_chi";
	int  L          = 15;
	double z_beg	= 40; 
	double z_end	= 40;
	double z_step	= 10;
	
	int myrank = 0;
	int nprocs = 1;
	
	int znum_total = (int)round((z_end - z_beg) / z_step) + 1;
	int nbeta = L + 1;
	
	int qnum = QNUM;
	double* qvals = mkarray(0.0, QMAX, qnum); 

	struct mol_prms *prms = mol_prms_read(prm_path);
	struct saxs_form_factor_table *ff_table = default_ff_table(map_path);
	
//=================================================================
//=========================== receptor ============================
	struct mol_atom_group *pa1 = mol_read_pdb(rec_path);
	mol_atom_group_add_prms(pa1, prms);
	
	struct mol_vector3 coe;
	center_of_extrema(&coe, pa1);
	MOL_VEC_MULT_SCALAR(coe, coe, -1.0);
	mol_atom_group_translate(pa1, &coe);

	struct sxs_spf_full* A = atom_grp2spf(pa1, ff_table, qvals, qnum, L, 1);

//=================================================================
//=========================== ligand ==============================
	struct mol_atom_group *pa2 = mol_read_pdb(lig_path);
	mol_atom_group_add_prms(pa2, prms);
	
	struct mol_vector3 com;
	centroid(&com, pa2);
	MOL_VEC_MULT_SCALAR(com, com, -1.0);
	mol_atom_group_translate(pa2, &com);
	
	struct sxs_spf_full* B = atom_grp2spf(pa2, ff_table, qvals, qnum, L, 1);

//=================================================================
//========================= experiment ============================
	struct sxs_profile* exp_profile = sxs_profile_read(exp_path);
	
	double mean_radius = (A->rm * pa1->natoms + B->rm * pa2->natoms) / 
	                     (pa1->natoms + pa2->natoms);
	                     
	struct sxs_opt_params* params = sxs_opt_params_create(exp_profile, qvals, qnum, mean_radius);
	
	mol_prms_free(prms);
	sxs_profile_free(exp_profile);
	
//=================================================================
//=========================== scoring =============================  
	
	double* zvals = calloc(znum_total, sizeof(double));
	int znum = 0;
	for (double zval = z_beg; zval < (z_end + 0.001); zval += z_step) {     
		zvals[znum++] = zval;
	}
	

	double* saxs_score;
	int*    saxs_index; 
	int*      ft_index;
	double*    c1_list;
	double*    c2_list;
	
	int ft_num = 0;    
	int ft_id;
	
	FILE* euler_file = fopen(eul_path, "r");
	
	double line[5];
	double z;
	while(fscanf(euler_file, "%d %lf %lf %lf %lf %lf %lf", &ft_id, &z, 
	             &line[0], &line[1], &line[2], &line[3], &line[4]) != EOF)
	{
		for (int j = 0; j < znum; j++) {
			if (zvals[j] > z - 0.001 && zvals[j] < z + 0.001) {
				ft_num++;
			}
		}
	}
	
	rewind(euler_file);

	saxs_index = calloc(ft_num, sizeof(int));
	saxs_score = calloc(ft_num, sizeof(double));
	ft_index = calloc(ft_num, sizeof(int));
	c1_list = calloc(ft_num, sizeof(double));
	c2_list = calloc(ft_num, sizeof(double));
	
	double b_step = M_PI / L;
	double a_step = 2.0 * M_PI / (2*L+1);
	int tmp, tmp_id;
	
	// convert euler coordinates into indices
	struct sxs_euler euler;
	int counter = 0;
	
	while(fscanf(euler_file, "%d %lf %lf %lf %lf %lf %lf", &ft_id, 
	             &euler.z,  &euler.b1, &euler.g1, &euler.a2, 
	             &euler.b2, &euler.g2) != EOF) {

		for (int j = 0; j < znum; j++) {
			if (zvals[j] > euler.z - 0.001 && zvals[j] < euler.z + 0.001) {
			
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

				saxs_index[counter] = tmp_id;
				ft_index[counter] = ft_id;
				counter++;
			}
		}
	}
	
	fclose(euler_file);
	
	compute_saxs_scores(saxs_score, 
	                    c1_list,
	                    c2_list, 
	                    saxs_index, 
	                    ft_num, 
	                    A, B, params, 
	                    qvals, qnum, 
	                    zvals, znum, 
	                    L, 1);

	FILE* refff = fopen("../ref_chi", "w");
	for (int i = 0; i < ft_num; i++) {
		fprintf(refff, "%d\t%.3lf\t%.3lf\t%.3lf\n", ft_index[i], saxs_score[i], c1_list[i], c2_list[i]);
	}
	fclose(refff);

	FILE* ref = fopen(ref_path, "r");
	
	double score, c1, c2;
	int nread;
	for (int i = 0; i < ft_num; i++) {
		nread = fscanf(ref, "%d %lf %lf %lf", &ft_id, &score, &c1, &c2);
		ck_assert_int_eq(nread, 4);
		ck_assert_int_eq(ft_id, ft_index[i]);
		ck_assert_double_eq_tol(score, saxs_score[i], tol1);
		ck_assert_double_eq_tol(c1, c1_list[i], tol1);
		ck_assert_double_eq_tol(c2, c2_list[i], tol1);
	}

	fclose(ref);
	
	sxs_myfree(saxs_score);
	sxs_myfree(saxs_index);
	sxs_myfree(ft_index);
	sxs_myfree(c1_list);
	sxs_myfree(c2_list);
	sxs_opt_params_free(params);
}
END_TEST

START_TEST(minimize_score)
{
	char map_path[] = "pdb_formfactor_mapping_clean.prm";          // for finding atoms in param file
	char prm_path[] = "atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec";          // param file
	char pdb_path[] = "4g9s_native_dimer.pdb";
	char exp_path[]	= "ref_saxs";
	char ref_path[] = "ref_fitted_profile";
	int  L = 15;
	int  qnum = 50;
	double* qvals = mkarray(0.0, 0.5, qnum); 

	struct mol_prms *prms = mol_prms_read(prm_path);
	struct saxs_form_factor_table *ff_table = default_ff_table(map_path);

	struct mol_atom_group *pa = mol_read_pdb(pdb_path);
	mol_atom_group_add_prms(pa, prms);

	struct mol_vector3 coe;
	center_of_extrema(&coe, pa);
	MOL_VEC_MULT_SCALAR(coe, coe, -1.0);
	mol_atom_group_translate(pa, &coe);

	struct sxs_spf_full* A = atom_grp2spf(pa, ff_table, qvals, qnum, L, 1);
	
	struct sxs_profile* exp_profile = sxs_profile_read(exp_path);
	struct sxs_opt_params* params = sxs_opt_params_create(exp_profile, qvals, qnum, A->rm);
	
	struct sxs_profile* profile = sxs_profile_create(qvals, qnum, 1);
	spf2fitted_profile(profile, A, params);
	sxs_profile_write("../ref_fitted_profile", profile);

	struct sxs_profile* ref_profile = sxs_profile_read(ref_path);
	cmp_profiles(profile, ref_profile);

	mol_prms_free(prms);
	sxs_profile_free(exp_profile);
	sxs_profile_free(profile);
	sxs_profile_free(ref_profile);
	sxs_opt_params_free(params);
	mol_atom_group_free(pa);
	sxs_myfree(ff_table);
}
END_TEST


Suite *featurex_suite(void)
{
	Suite *suite = suite_create("featurex");


	// Add test cases here
	// Each test case can call multiple test functions
	// Too add a test case, call tcase_add_test
	// The first argument is the TCase struct, the second is the
	//  test function name.
	TCase *tcase1 = tcase_create("test_pdb2spf");
	tcase_add_test(tcase1, test_pdb2spf);
	
	TCase *tcase2 = tcase_create("test_sxs_profile_from_spf");
	tcase_add_test(tcase2, test_sxs_profile_from_spf);
	
	TCase *tcase3 = tcase_create("test_main_fun");
	tcase_set_timeout(tcase3, 100);
	tcase_add_test(tcase3, score_conformations);
	
	TCase *tcase4 = tcase_create("test_minimizer");
	tcase_add_test(tcase4, minimize_score);
	
	suite_add_tcase(suite, tcase1);
	suite_add_tcase(suite, tcase2);
	suite_add_tcase(suite, tcase3);
	suite_add_tcase(suite, tcase4);

	return suite;
}

int main(void)
{
	Suite *suite = featurex_suite();
	SRunner *runner = srunner_create(suite);
	srunner_run_all(runner, CK_NORMAL);

	int number_failed = srunner_ntests_failed(runner);
	srunner_free(runner);
	return number_failed;
}

