#ifndef _FORM_FACTOR_H_
#define _FORM_FACTOR_H_

#include "mol2/atom_group.h"
#include "saxs_utils.h"

#include "khash.h"


struct saxs_resi_atom {
	char residue_name[8];
	char atom_name[8];
};

KHASH_DECLARE(SAXS_FF_HASH, struct saxs_resi_atom, enum saxs_ff_type);

enum saxs_ff_type
{
	H, He,
	Li, Be, B, C, N, O, F, Ne, // periodic table, lines 1-2 (10)
	Na, Mg, Al, Si, P, S, Cl, Ar, // line 3 (8)
	K, Ca, Cr, Mn, Fe, Co,
	Ni, Cu, Zn, Se, Br, // line 4 (11)
	Io, Ir, Pt, Au, Hg, SINGLE_ATOM_SIZE = 34,
	CH=34, CH2=35, CH3=36, NH=37, NH2=38, NH3=39, OH=40, s_OH2=41, SH=42, PO4=43,
	HEAVY_ATOM_SIZE=44, s_UNK=99
};

// f[q] = c + (sum_i a_i*e^(-b_i*q(q^2)))
struct saxs_form_factor_coeffs
{
	double a[5];
	double b[5];
	double c;
	double excl_vol;

	bool valid;
};

struct saxs_form_factor
{
	// Need vacuum_ff, dummy_ff
	double zero_ff;
	double vacuum_ff;
	double dummy_ff;

	double *values;
};

struct saxs_form_factor_table
{
	struct saxs_form_factor factors[HEAVY_ATOM_SIZE];
	struct saxs_form_factor_coeffs coeffs[HEAVY_ATOM_SIZE];

	khash_t(SAXS_FF_HASH) *type_map;

	// Temp variables
	struct saxs_form_factor *water_ff;
};


struct saxs_form_factor_table *default_ff_table(const char *type_mapping_file);
struct saxs_form_factor_table *read_ff_table(const char *ff_file, const char *type_mapping_file);
void ff_table_destroy(struct saxs_form_factor_table *table);

const char *ff_type_to_string(enum saxs_ff_type type);

const struct saxs_form_factor *get_ff(const struct saxs_form_factor_table *table,
				      const struct mol_atom_group *ag, size_t atom_index);

double dummy_ff(const struct saxs_form_factor_table *table,
		const struct mol_atom_group *ag, size_t atom_index);
double vacuum_ff(const struct saxs_form_factor_table *table,
		 const struct mol_atom_group *ag, size_t atom_index);
double compute_form_factor(const struct saxs_form_factor_table *table, enum saxs_ff_type type, double q);
void build_form_factor_tables(struct saxs_form_factor_table *table,
                              double *q, size_t n_q);
enum sax_ff_type atom2ff_type(const struct saxs_form_factor_table *table,
			      const struct mol_atom_group *ag, size_t atom_index);

void get_dummy_ff(double *dummy_ffs,
                  const struct mol_atom_group *ag,
                  const struct saxs_form_factor_table *table);
void get_vacuum_ff(double *vacuum_ffs,
                   const struct mol_atom_group *ag,
		   const struct saxs_form_factor_table *table);

// debugging functions
void dump_mappings(struct saxs_form_factor_table *table);
void print_ff_table(struct saxs_form_factor_table *table);

#endif /* _FORM_FACTOR_H_ */
