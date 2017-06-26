/** @file
 * \brief Contains functions for convertation of ft and rm files
 * to Euler angles.
 * Center of mass - COM.
 * Center of extrema - COE.
 */
 
/** \addtogroup sph_interface
 * @{
 */
#pragma once

#include "common.h"
#include "saxs_utils.h"

#include "mol2/vector.h"
#include "mol2/matrix.h"
#include "mol2/lists.h"

/** 
 * Contains euler coordinates of a conformation.
 * z - translation of the COM of the ligand with 
 * respect to the COE of the receptor.
 */
struct sxs_euler
{
	double z;
	double b1;
	double g1;
	double a2;
	double b2;
	double g2;
};

struct sxs_index 
{
	int z ;
	int b1;
	int g1;
	int a2;
	int b2;
	int g2;
};

size_t sxs_assemble_index (struct sxs_index* id, int nbeta, int L);

void sxs_disassemble_index (struct sxs_index* id, size_t index, int nbeta, int L);

/**
 * Convert dimer conformation represented by the translation of the ligand's COM and 
 * rotation matrix around its center.
 *
 * @param[out] euler Euler coordinates.
 * @param tv Translation of ligand's COM with respect to the receptor's COM.
 * @param rm Rotation of the ligand around its COM.
 * @param ref_lig Reference position of the ligand's COM with respect to the receptor's COE
 * in the pdb file.
 */
void sxs_ft2euler (struct sxs_euler* euler, struct mol_vector3 *tv, struct mol_matrix3 *rm, struct mol_vector3 *ref_lig);

/**
 * Convert ft file and rm file (file with rotation matrices) to Euler coordinates file.
 * @param eu_path Output path.
 * @param ft_path FT-file.
 * @param rm_path Rotation matrices file.
 * @param ref_lig Reference position of the ligand's COM with respect to the receptor's COE.
 */
void sxs_ft_file2euler_file (const char* eu_path, const char* ft_path, const char* rm_path, struct mol_vector3 *ref_lig);

/** @} */


