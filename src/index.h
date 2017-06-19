#pragma once

#include "common.h"

#include "mol2/vector.h"
#include "mol2/matrix.h"
#include "mol2/lists.h"

struct saxs_euler
{
	double z;
	double b1;
	double g1;
	double a2;
	double b2;
	double g2;
};

struct saxs_index 
{
	int z ;
	int b1;
	int g1;
	int a2;
	int b2;
	int g2;
};

size_t saxs_assemble_index (struct saxs_index* id, int nbeta, int L);

void saxs_disassemble_index (struct saxs_index* id, size_t index, int nbeta, int L);

void saxs_fill_active_rotation_matrix(struct mol_matrix3 *rm, double alpha, double beta, double gamma);

void saxs_mult_rot_mats(struct mol_matrix3 *c, struct mol_matrix3 *a, struct mol_matrix3 *b);

void saxs_ft2euler (struct saxs_euler* euler, struct mol_vector3 *tv, struct mol_matrix3 *rm, struct mol_vector3 *ref_lig);

void saxs_ft_file2euler_file (const char* eu_path, const char* ft_path, const char* rm_path, struct mol_vector3 *ref_lig);




