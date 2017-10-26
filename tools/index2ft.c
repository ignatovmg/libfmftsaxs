
#include "common.h"

#include "mol2/pdb.h"
#include "mol2/atom_group.h"
#include "mol2/vector.h"
#include "mol2/matrix.h"

#include "index.h"

int main(int argc, char** argv)
{
	int opt = 1;
	char* rec_path = argv[opt++];
	char* lig_path = argv[opt++];
	char* eu_path  = argv[opt++];
	char* ft_path  = argv[opt++];
	char* rm_path  = argv[opt++];

	struct mol_atom_group* rec = mol_read_pdb(rec_path);
	struct mol_atom_group* lig = mol_read_pdb(lig_path);
	
	struct mol_vector3 coe;
	center_of_extrema(&coe, rec);

	struct mol_vector3 com;
	centroid(&com, lig);
	
	struct mol_vector3 ref_lig;
	MOL_VEC_SUB(ref_lig, com, coe);

	mol_atom_group_free(rec);
	mol_atom_group_free(lig);
	
	FILE* eu_file = fopen(eu_path, "r");
	FILE* ft_file = fopen(ft_path, "w");
	FILE* rm_file = fopen(rm_path, "w");
	
	int L;
	double z, b1, b2;
	if (fscanf(eu_file, "%d %lf %lf %lf", &L, &z, &b1, &b2) != 4) {
		ERROR_MSG("Wrong input file format");
	}
	
	int iline = 0;
	double score, c1, c2;
	double astep = 2*M_PI / (2*L+1);
	struct mol_matrix3 pm, am, rm;
	struct mol_vector3 tr;

	for (int ia2 = 0; ia2 < 2*L+1; ia2++) {
		for (int ig1 = 0; ig1 < 2*L+1; ig1++) {
			for (int ig2 = 0; ig2 < 2*L+1; ig2++) {
			
				if (fscanf(eu_file, "%lf %lf %lf", &score, &c1, &c2) != 3) {
					ERROR_MSG("Wrong input file format");
				}
				
				sxs_fill_passive_rotation_matrix(&pm, 0, b1, ig1*astep);			
				tr.X = pm.m13 * z;
				tr.Y = pm.m23 * z;
				tr.Z = pm.m33 * z;

				sxs_fill_active_rotation_matrix(&am, -ia2*astep, b2, -ig2*astep);
				sxs_mult_rot_mats(&rm, &pm, &am);

				fprintf(rm_file, "%d \t % .9lf % .9lf % .9lf % .9lf % .9lf % .9lf % .9lf % .9lf % .9lf\n",
				        iline, rm.m11, rm.m12, rm.m13, rm.m21, rm.m22, 
				        rm.m23, rm.m31, rm.m32, rm.m33);
				
				fprintf(ft_file, "%d \t %.2lf %.2lf %.2lf %.3lf 0.0 0.0 0.0 0.0 0.0\n", 
				        iline, tr.X - ref_lig.X, tr.Y - ref_lig.Y, tr.Z - ref_lig.Z, score);

				iline++;
			}
		}
	}

	fclose(eu_file);
	fclose(ft_file);
	fclose(rm_file);


	return EXIT_SUCCESS;
}
