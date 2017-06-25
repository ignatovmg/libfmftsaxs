
#include "common.h"

#include "mol2/pdb.h"
#include "mol2/atom_group.h"
#include "mol2/vector.h"

#include "index.h"

int main(int argc, char** argv)
{
	char* ft_path = argv[1];
	char* rm_path = argv[2];
	char* rec_path = argv[3];
	char* lig_path = argv[4];
	char* eu_path = argv[5];

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

	sxs_ft_file2euler_file (eu_path, ft_path, rm_path, &ref_lig);

	return EXIT_SUCCESS;
}
