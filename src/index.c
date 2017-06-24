#include "index.h"

size_t sxs_assemble_index (struct sxs_index* id, int nbeta, int L)
{
	int fft = 2*L + 1;
	return ((((id->z * nbeta + id->b1) * nbeta + id->b2) * fft + id->a2) * fft + id->g1) * fft + id->g2;
}

void sxs_disassemble_index (struct sxs_index* id, size_t index, int nbeta, int L)
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

static double check_range(double a)
{
	if (a < -1.0) { return -1.0; }
	if (a >  1.0) { return  1.0; }
	return a;
}

void sxs_ft2euler (struct sxs_euler* euler, struct mol_vector3 *tv, struct mol_matrix3 *rm, struct mol_vector3 *ref_lig)
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
	sxs_fill_active_rotation_matrix(&rec_rm, 0.0, b1, g1);
	sxs_mult_rot_mats(&lig_rm, &rec_rm, rm);
	
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

void sxs_ft_file2euler_file (const char* eu_path, const char* ft_path, const char* rm_path, struct mol_vector3 *ref_lig)
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
	
	struct sxs_euler euler;
	struct mol_vector3 vec;
	struct mol_matrix3_list* mol_mat_list = mol_matrix3_list_from_file(rm_path);
	struct mol_matrix3* mat_list = mol_mat_list->members;
	
	int nread;
	for (int i = 0; i < ln; i++) {
		nread = fscanf(ft_file, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", \
		       &ft_id, &ft_line[0], &ft_line[1], &ft_line[2], &ft_line[3], \
		       &ft_line[4], &ft_line[5], &ft_line[6], &ft_line[7], &ft_line[8]);
		       
		if (nread != 10) {
			ERROR_MSG("Wrong input file format.");
		}
	
		vec.X = ft_line[0];
		vec.Y = ft_line[1];
		vec.Z = ft_line[2];
		
		sxs_ft2euler(&euler, &vec, &mat_list[ft_id], ref_lig);
		
		fprintf(eu_file, "%d\t% .3f\t% .3f\t% .3f\t% .3f\t% .3f\t% .3f\n", 
		       ft_id, euler.z, euler.b1, euler.g1,
		       euler.a2, euler.b2, euler.g2);
	}
	
	fclose(eu_file);
	fclose(ft_file);
}




