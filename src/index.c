#include "index.h"

size_t saxs_assemble_index (struct saxs_index* id, int nbeta, int L)
{
	int fft = 2*L + 1;
	return ((((id->z * nbeta + id->b1) * nbeta + id->b2) * fft + id->a2) * fft + id->g1) * fft + id->g2;
}

void saxs_disassemble_index (struct saxs_index* id, size_t index, int nbeta, int L)
{
	int fft = 2 * L + 1;
	int tmp;
	
	id->g2 = index % fft;
	tmp = index / fft;
	
	id->g1 = tmp % fft;
	tmp = tmp / fft;
	
	id->a2 = tmp % fft;
	tmp = tmp / fft;
	
	id->b2 = tmp % nbeta;
	tmp = tmp / nbeta;
	
	id->b1 = tmp % nbeta;
	
	id->z  = tmp / nbeta;
}


void saxs_fill_active_rotation_matrix(struct mol_matrix3 *rm, double alpha, double beta, double gamma)
{
	rm->m11 = cos(gamma)*cos(beta)*cos(alpha) - sin(gamma)*sin(alpha);
	rm->m21 = cos(gamma)*cos(beta)*sin(alpha) + sin(gamma)*cos(alpha);
	rm->m31 = -cos(gamma)*sin(beta);
	rm->m12 = -sin(gamma)*cos(beta)*cos(alpha) - cos(gamma)*sin(alpha);
	rm->m22 = -sin(gamma)*cos(beta)*sin(alpha) + cos(gamma)*cos(alpha);
	rm->m32 = sin(gamma)*sin(beta);
	rm->m13 = sin(beta)*cos(alpha);
	rm->m23 = sin(beta)*sin(alpha);
	rm->m33 = cos(beta);
//	{ cos(g)*cos(b)*cos(a) - sin(g)*sin(a),  cos(g)*cos(b)*sin(a) + sin(g)*cos(a), -cos(g)*sin(b)},
//	{-sin(g)*cos(b)*cos(a) - cos(g)*sin(a), -sin(g)*cos(b)*sin(a) + cos(g)*cos(a),  sin(g)*sin(b)},
//	{ sin(b)*cos(a)                       ,  sin(b)*sin(a)                       ,  cos(b)       }
}

void saxs_mult_rot_mats(struct mol_matrix3 *c, struct mol_matrix3 *a, struct mol_matrix3 *b)
{
	c->m11 = a->m11 * b->m11 + a->m12 * b->m21 + a->m13 * b->m31;
	c->m12 = a->m11 * b->m12 + a->m12 * b->m22 + a->m13 * b->m32;
	c->m13 = a->m11 * b->m13 + a->m12 * b->m23 + a->m13 * b->m33;

	c->m21 = a->m21 * b->m11 + a->m22 * b->m21 + a->m23 * b->m31;
	c->m22 = a->m21 * b->m12 + a->m22 * b->m22 + a->m23 * b->m32;
	c->m23 = a->m21 * b->m13 + a->m22 * b->m23 + a->m23 * b->m33;

	c->m31 = a->m31 * b->m11 + a->m32 * b->m21 + a->m33 * b->m31;
	c->m32 = a->m31 * b->m12 + a->m32 * b->m22 + a->m33 * b->m32;
	c->m33 = a->m31 * b->m13 + a->m32 * b->m23 + a->m33 * b->m33;
}

static double check_range(double a)
{
	if (a < -1.0) { return -1.0; }
	if (a >  1.0) { return  1.0; }
	return a;
}

void saxs_ft2euler (struct saxs_euler* euler, struct mol_vector3 *tv, struct mol_matrix3 *rm, struct mol_vector3 *ref_lig)
{	
	struct mol_vector3 vec;
	MOL_VEC_ADD(vec, *tv, *ref_lig);
	
	double z  = round(sqrt(MOL_VEC_SQ_NORM(vec)));
	double b1 = acos(check_range(vec.Z / z));
	double g1 = acos(check_range( -vec.X / (z * sin(b1))));
	
	if (vec.Y / (z * sin(b1)) < 0.0) {
		g1 = 2*M_PI - g1;
	}
		
	struct mol_matrix3 rec_rm;
	struct mol_matrix3 lig_rm;
	saxs_fill_active_rotation_matrix(&rec_rm, 0.0, b1, g1);
	saxs_mult_rot_mats(&lig_rm, &rec_rm, rm);
	
	double b2 = acos(check_range(lig_rm.m33));
	double a2 = acos(check_range(lig_rm.m13 / sin(b2)));
	
	if (lig_rm.m23 / sin(b2) < 0.0) {
		a2 = 2*M_PI - a2;
	}
		
	double g2 = acos(check_range(-lig_rm.m31 / sin(b2)));
	
	if (lig_rm.m32 / sin(b2) < 0.0) {
		g2 = 2*M_PI - g2;
	}
	
	euler->z  = z;
	euler->b1 = b1;
	euler->g1 = g1;
	euler->a2 = a2;
	euler->b2 = b2;
	euler->g2 = g2;
}

void saxs_ft_file2euler_file (const char* eu_path, const char* ft_path, const char* rm_path, struct mol_vector3 *ref_lig)
{
	int    ft_id, ln;
	double ft_line[9];
	
	ln = 0;
	
	FILE* ft_file = fopen(ft_path, "r");
	FILE* eu_file = fopen(eu_path, "w");
	
	while(fscanf(ft_file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", \
		&ft_id, &ft_line[0], &ft_line[1], &ft_line[2], &ft_line[3], &ft_line[4], \
		&ft_line[5], &ft_line[6], &ft_line[7], &ft_line[8]) != EOF) { ln++; }
		
	rewind(ft_file);
	
	//struct mol_vector3 *vec_list = calloc(ln, sizeof(struct mol_vector3));
	//int *ft_id_list = calloc(ln, sizeof(int));
	//int *saxs_id_list = calloc(ln, sizeof(int));
	
	struct saxs_euler euler;
	struct mol_vector3 vec;
	struct mol_matrix3_list* mol_mat_list = mol_matrix3_list_from_file(rm_path);
	struct mol_matrix3* mat_list = mol_mat_list->members;
	
	for (int i = 0; i < ln; i++) {
		fscanf(ft_file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", \
		       &ft_id, &ft_line[0], &ft_line[1], &ft_line[2], &ft_line[3], \
		       &ft_line[4], &ft_line[5], &ft_line[6], &ft_line[7], &ft_line[8]);
	
		vec.X = ft_line[0];
		vec.Y = ft_line[1];
		vec.Z = ft_line[2];
		
		saxs_ft2euler(&euler, &vec, &mat_list[ft_id], ref_lig);
		
		fprintf(eu_file, "%d\t% .3f\t% .3f\t% .3f\t% .3f\t% .3f\t% .3f\n", 
		       ft_id, euler.z, euler.b1, euler.g1,
		       euler.a2, euler.b2, euler.g2);
	}
	
	fclose(eu_file);
	fclose(ft_file);
	
	/*int tmp_id, tmp;
	double beta_step = M_PI / (nbeta - 1);
	double alpha_step = 2*M_PI / (2*L+1);
	struct mol_vector3* mat_list = mol_matrix3_list_from_file(rm_path);*/
	
	
	/*for (int i = 0; i < ln; i++) {
	
		struct saxs_euler euler;
		
		saxs_ft2euler(&euler, &vec_list[i], &mat_list[i], ref_lig);

		if (euler->z < z_beg || euler->z > z_end) {
			continue;
		}
		
		//
		euler->a2 = 2*M_PI - euler->a2;
		euler->g2 = 2*M_PI - euler->g2;
	
        tmp_id = (int)(round((euler->z - z_beg) / z_step)) * nbeta;

        tmp = (int)(round(euler->b1 / beta_step));
        tmp_id = (tmp_id + tmp) * fragnum;

        tmp = (int)(round(euler->b2 / beta_step));
        tmp_id = (tmp_id + tmp) * (2*L+1);

        tmp = (int)(round(euler->a2 / alpha_step));
        tmp_id = (tmp_id + tmp) * (2*L+1);

        tmp = (int)(round(euler->g1 / alpha_step));
        tmp_id = (tmp_id + tmp) * (2*L+1);

        tmp = (int)(round(euler->g2 / alpha_step));
        tmp_id = tmp_id + tmp;

        saxs_id_list[i] = tmp_id;
	}*/
}




