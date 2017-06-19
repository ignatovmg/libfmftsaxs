#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#include "form_factor_table.h"

#include "saxs_utils.h"

static double default_zero_form_factors[] = {
	-0.720147, -0.720228,
	// H          He
	1.591, 2.591, 3.591, 0.50824, 6.16294, 4.94998, 7.591, 6.993,
	// Li    Be     B      C        N        O        F      Ne
	7.9864, 8.9805, 9.984, 10.984, 13.0855, 9.36656, 13.984, 16.591,
	// Na     Mg      Al      Si      P       S         Cl      Ar
	15.984, 14.9965, 20.984, 21.984, 20.9946, 23.984,
	// K       Ca2+     Cr      Mn      Fe2+     Co
	24.984, 25.984, 24.9936, 30.9825, 31.984, 49.16,
	// Ni      Cu      Zn2+     Se       Br      I
	70.35676, 71.35676, 72.324, 73.35676,
	// Ir        Pt        Au      Hg
	-0.211907, -0.932054, -1.6522, 5.44279, 4.72265, 4.0025, 4.22983, 3.50968, 8.64641
	// CH         CH2        CH3     NH       NH2      NH3     OH       OH2      SH
};

static double default_vacuum_zero_form_factors[] = {
	0.999953, 0.999872, 2.99, 3.99, 4.99, 5.9992, 6.9946, 7.9994, 8.99, 9.999,
	// H        He        Li    Be    B     C       N       O       F     Ne
	10.9924, 11.9865, 12.99, 13.99, 14.9993, 15.9998, 16.99, 17.99,
	// Na       Mg       Al     Si     P        S        Cl     Ar
	18.99, 18.0025,  23.99, 24.99,  24.0006, 26.99,
	// K      Ca2+      Cr     Mn      Fe2+     Co
	27.99, 28.99, 27.9996, 33.99, 34.99, 52.99, 76.99, 77.99, 78.9572, 79.99,
	// Ni     Cu     Zn2+     Se     Br     I      Ir     Pt     Au       Hg
	6.99915, 7.99911, 8.99906, 7.99455, 8.99451, 9.99446, 8.99935, 9.9993, 16.9998
	// CH      CH2      CH3      NH       NH2      NH3      OH       OH2      SH
};

static double default_dummy_zero_form_factors[] = {
	1.7201, 1.7201, 1.399, 1.399, 1.399 , 5.49096, 0.83166, 3.04942, 1.399, 3.006,
	// H      He      Li?    Be?    B?      C        N        O        F?     Ne
	3.006, 3.006, 3.006, 3.006, 1.91382, 6.63324, 3.006, 1.399,
	// Na    Mg     Al?    Si?    P        S        Cl?    Ar?
	3.006, 3.006, 3.006, 3.006, 3.006, 3.006,
	// K?    Ca2+   Cr?    Mn?    Fe2+   Co?
	3.006, 3.006, 3.006, 3.006, 3.006, 3.83, 6.63324, 6.63324, 6.63324, 6.63324,
	// Ni?   Cu?    Zn2+   Se     Br?    I?    Ir?      Pt?      Au       Hg
	7.21106, 8.93116, 10.6513, 2.55176, 4.27186, 5.99196, 4.76952, 6.48962, 8.35334
	// CH      CH2       CH3     NH       NH2      NH3      OH       OH2      SH
};

// Taken from khash.h
static kh_inline khint_t __saxs_resi_atom_hash_func(const struct saxs_resi_atom key)
{
	khint_t h = (khint_t) key.residue_name[0];
	for (size_t i = 1; key.residue_name[i]; i++) {
		h = (h << 5) - h + (khint_t) key.residue_name[i];
	}
	for (size_t i = 0; key.atom_name[i]; i++) {
		h = (h << 5) - h + (khint_t) key.atom_name[i];
	}

	return h;
}

static khint_t __saxs_resi_atom_equal(const struct saxs_resi_atom a,
				      const struct saxs_resi_atom b)
{
	return (strncmp(a.residue_name, b.residue_name, 7) == 0 &&
		strncmp(a.atom_name, b.atom_name, 7) == 0);
}

__KHASH_IMPL(SAXS_FF_HASH, , struct saxs_resi_atom, enum saxs_ff_type, 1,
	     __saxs_resi_atom_hash_func, __saxs_resi_atom_equal);


static enum saxs_ff_type string_to_ff_type(const char *name)
{
	static char cmp[8];
	strcpy(cmp, name);
	upcase(cmp);

	if (!strcmp(name, "H"))     return H;
	if (!strcmp(name, "HE"))    return He;
	if (!strcmp(name, "C"))     return C;
	if (!strcmp(name, "N"))     return N;
	if (!strcmp(name, "O"))     return O;
	if (!strcmp(name, "NE"))    return Ne;
	if (!strcmp(name, "SOD+"))  return Na;
	if (!strcmp(name, "MG2+"))  return Mg;
	if (!strcmp(name, "P"))     return P;
	if (!strcmp(name, "S"))     return S;
	if (!strcmp(name, "K"))     return K;
	if (!strcmp(name, "CAL2+")) return Ca;
	if (!strcmp(name, "FE2+"))  return Fe;
	if (!strcmp(name, "ZN2+"))  return Zn;
	if (!strcmp(name, "SE"))    return Se;
	if (!strcmp(name, "AU"))    return Au;
	if (!strcmp(name, "CH"))    return CH;
	if (!strcmp(name, "CH2"))   return CH2;
	if (!strcmp(name, "CH3"))   return CH3;
	if (!strcmp(name, "NH"))    return NH;
	if (!strcmp(name, "NH2"))   return NH2;
	if (!strcmp(name, "NH3"))   return NH3;
	if (!strcmp(name, "OH"))    return OH;
	if (!strcmp(name, "SH"))    return SH;
	return s_UNK;
}

const char *ff_type_to_string(enum saxs_ff_type type)
{
	if (type == H) return "H";
	if (type == He) return "HE";
	if (type == C) return "C";
	if (type == N) return "N";
	if (type == O) return "O";
	if (type == Ne) return "NE";
	if (type == Na) return "SOD+";
	if (type == Mg) return "MG2+";
	if (type == P) return "P";
	if (type == S) return "S";
	if (type == K) return "K";
	if (type == Ca) return "CAL2+";
	if (type == Fe) return "FE2+";
	if (type == Zn) return "ZN2+";
	if (type == Se) return "SE";
	if (type == Au) return "AU";
	if (type == CH) return "CH";
	if (type == CH2) return "CH2";
	if (type == CH3) return "CH3";
	if (type == NH) return "NH";
	if (type == NH2) return "NH2";
	if (type == NH3) return "NH3";
	if (type == OH) return "OH";
	if (type == SH) return "SH";
	return NULL;
}

struct saxs_form_factor_table *default_ff_table(const char *type_mapping_file)
{
	static struct saxs_form_factor_table *table = NULL;

	if (table != NULL)
		return table;

	FILE *tm_fp = fopen(type_mapping_file, "r");
	if (tm_fp == NULL)
	{
		fprintf(stderr, "Could not open file %s\n", type_mapping_file);
		exit(1);
	}

	table = calloc(1, sizeof(struct saxs_form_factor_table));
	table->type_map = kh_init(SAXS_FF_HASH);
	for (int i = 0; i < HEAVY_ATOM_SIZE; ++i) {
		table->factors[i].zero_ff = default_zero_form_factors[i];
		table->factors[i].vacuum_ff = default_vacuum_zero_form_factors[i];
		table->factors[i].dummy_ff = default_dummy_zero_form_factors[i];
	}

	char *line = NULL;
	size_t len;

	struct saxs_resi_atom line_key = {{0}, {0}};
	char ff_type_name[16];
	enum saxs_ff_type ff_type;
	int ret;

	while (getline(&line, &len, tm_fp) != -1)
	{
		if (line[0] == '#' || is_whitespace_line(line)) {
			continue;
		}

		//            res  atom  ff_type
		sscanf(line, "%7s  %7s   %15s", line_key.residue_name, line_key.atom_name, ff_type_name);
		ff_type = string_to_ff_type(ff_type_name);

		khiter_t val = kh_put(SAXS_FF_HASH, table->type_map, line_key, &ret);
		if (val != kh_end(table->type_map)) {
			kh_value(table->type_map, val) = ff_type;
		}
	}
	fclose(tm_fp);

	free_if_not_null(line);

	return table;
}

struct saxs_form_factor_table *read_ff_table(const char *ff_file, const char *type_mapping_file)
{
	FILE *ff_fp = fopen(ff_file, "r");
	if (ff_fp == NULL) {
		return NULL;
	}

	struct saxs_form_factor_table *table = calloc(1, sizeof(struct saxs_form_factor_table));

	char *line = NULL;
	size_t len;

	char group_name[16];
	double a1, a2, a3, a4, a5;
	double b1, b2, b3, b4, b5;
	double c, excl_vol;

	struct saxs_form_factor_coeffs *cc;

	while (getline(&line, &len, ff_fp) != -1)
	{
		if (line[0] == '#' || is_whitespace_line(line)) {
			continue;
		}

		sscanf(line,
		       //gn   a1  a2  a3  a4  a5  c   b1  b2  b3  b4  b5  vol
		       "%15s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		       group_name,
		       &a1, &a2, &a3, &a4, &a5,
		       &c,
		       &b1, &b2, &b3, &b4, &b5,
		       &excl_vol);
		enum saxs_ff_type type = string_to_ff_type(group_name);
		if (type != s_UNK)
		{
			// Fill in coefficient data
			cc = &(table->coeffs[type]);
			cc->a[0] = a1; cc->a[1] = a2; cc->a[2] = a3; cc->a[3] = a4; cc->a[4] = a5;
			cc->b[0] = b1; cc->b[1] = b2; cc->b[2] = b3; cc->b[3] = b4; cc->b[4] = b5;
			cc->excl_vol = excl_vol;
			cc->valid = true;
		}
	}
	fclose(ff_fp);

	// Read type mappings
	FILE *tm_fp = fopen(type_mapping_file, "r");

	struct saxs_resi_atom line_key = {{0}, {0}};
	char ff_type_name[16];
	enum saxs_ff_type ff_type;
	int ret;

	while (getline(&line, &len, tm_fp) != -1)
	{
		if (line[0] == '#' || is_whitespace_line(line)) {
			continue;
		}

		//            res  atom  ff_type
		sscanf(line, "%7s  %7s   %15s", line_key.residue_name, line_key.atom_name, ff_type_name);
		ff_type = string_to_ff_type(ff_type_name);

		khiter_t val = kh_put(SAXS_FF_HASH, table->type_map, line_key, &ret);
		if (val != kh_end(table->type_map)) {
			kh_value(table->type_map, val) = ff_type;
		}
	}
	fclose(tm_fp);

	free_if_not_null(line);

	return table;
}

void ff_table_destroy(struct saxs_form_factor_table *table)
{
	if (table != NULL)
	{
		if (table->coeffs != NULL)
		{
			for (int i = 0; i < HEAVY_ATOM_SIZE; ++i)
			{
				struct saxs_form_factor_coeffs *coeffs = &(table->coeffs[i]);
				free_if_not_null(coeffs->a);
				free_if_not_null(coeffs->b);
			}
			free(table->coeffs);
		}

		if (table->factors != NULL)
		{
			for (int i = 0; i < HEAVY_ATOM_SIZE; ++i)
			{

			}
			free(table->factors);
		}

		free(table);
	}
}

const struct saxs_form_factor *get_ff(const struct saxs_form_factor_table *table,
				      const struct mol_atom_group *ag, size_t atom_index)
{
	struct saxs_resi_atom search_key = {{0}, {0}};
	char *src, *dst;
	src = ag->residue_name[atom_index];
	dst = search_key.residue_name;

	while(*src && *src == ' ') {
		src++;
	}
	while(*src && *src != ' ') {
		*dst = *src;
		dst++;
		src++;
	}

	src = ag->atom_name[atom_index];
	dst = search_key.atom_name;
	while(*src && *src == ' ') {
		src++;
	}
	while(*src && *src != ' ') {
		*dst = *src;
		dst++;
		src++;
	}

	khiter_t key;
	key = kh_get(SAXS_FF_HASH, table->type_map, search_key);

	if (key == kh_end(table->type_map)) {
		fprintf(stderr, "mapping for (%s, %s) not found\n",
			search_key.residue_name, search_key.atom_name);
		return NULL;
	} else {
		enum saxs_ff_type ff_type = kh_value(table->type_map, key);
		if (ff_type != s_UNK) {
			return &(table->factors[ff_type]);
		} else {
			printf("(%s, %s) maps to UNK\n",
			       search_key.residue_name, search_key.atom_name);
			return NULL;
		}
	}
}

double compute_form_factor(const struct saxs_form_factor_table *table, enum saxs_ff_type type, double q)
{
	const struct saxs_form_factor_coeffs *coeff = &(table->coeffs[type]);

	double f = coeff->c;
	// f = c + (SUM_i a_i*e^(-b_i*(q^2)))
	for (int i = 0; i < 5; ++i) {
		f += coeff->a[i] * exp(-coeff->b[i] * (q*q));
	}

	return f;
}

void build_form_factor_tables(struct saxs_form_factor_table *table, double *q_vals,
                              size_t n_qvals)
{
	static double two_thirds = 2.0/3.0;
	static double one_over_four_pi = 1.0/(4.0*M_PI);
	static double rho = 0.334; // TODO: replace with electron density of solvent

	// Precompute qq and ss;
	double *qq = calloc(n_qvals, sizeof(double));
	for (size_t i = 0; i < n_qvals; i++) {
		qq[i] = q_vals[i] * q_vals[i];
	}
	double *ss = (double *) malloc(n_qvals * sizeof(double));
	for (size_t i = 0; i < n_qvals; i++) {
		ss[i] = qq[i] * one_over_four_pi * one_over_four_pi;
	}

	// Create tables for single atoms
	for (int ia = 0; ia < SINGLE_ATOM_SIZE; ia++) {
		struct saxs_form_factor_coeffs *coeff = &(table->coeffs[ia]);
		if (!coeff->valid)
			continue;
		struct saxs_form_factor *factor = &(table->factors[ia]);

		for (size_t iq = 0; iq < n_qvals; iq++) {
			factor->values[iq] = coeff->c;
			// f = c + (SUM_i a_i*e^(-b_i*(q^2)))
			for (int j = 0; j < 5; j++) {
				factor->values[iq] += coeff->a[j] * exp(-coeff->b[j] * (ss[iq]));
			}

			double vol_c = - pow(coeff->excl_vol, two_thirds) * one_over_four_pi;
			factor->values[iq] -= rho * coeff->excl_vol * exp(vol_c * qq[iq]);
		}
		factor->zero_ff = coeff->c;
		for (int j = 0; j < 5; j++) {
			factor->zero_ff += coeff->a[j];
		}
		factor->vacuum_ff = factor->zero_ff;
		factor->dummy_ff = rho * coeff->excl_vol;
		factor->zero_ff -= rho * coeff->excl_vol;
	}

	// Create tables for compound groups
	struct saxs_form_factor *CH_factor = &(table->factors[CH]);
	struct saxs_form_factor *CH2_factor = &(table->factors[CH2]);
	struct saxs_form_factor *CH3_factor = &(table->factors[CH3]);
	struct saxs_form_factor *NH_factor = &(table->factors[NH]);
	struct saxs_form_factor *NH2_factor = &(table->factors[NH2]);
	struct saxs_form_factor *NH3_factor = &(table->factors[NH3]);
	struct saxs_form_factor *OH_factor = &(table->factors[OH]);
	struct saxs_form_factor *OH2_factor = &(table->factors[s_OH2]);
	struct saxs_form_factor *SH_factor = &(table->factors[SH]);
	struct saxs_form_factor *PO4_factor = &(table->factors[PO4]);

	CH_factor->values = calloc(n_qvals, sizeof(double));
	CH2_factor->values = calloc(n_qvals, sizeof(double));
	CH3_factor->values = calloc(n_qvals, sizeof(double));
	NH_factor->values = calloc(n_qvals, sizeof(double));
	NH2_factor->values = calloc(n_qvals, sizeof(double));
	NH3_factor->values = calloc(n_qvals, sizeof(double));
	OH_factor->values =  calloc(n_qvals, sizeof(double));
	OH2_factor->values = calloc(n_qvals, sizeof(double));
	SH_factor->values = calloc(n_qvals, sizeof(double));
	PO4_factor->values = calloc(n_qvals, sizeof(double));

	// Fill in factor table by summing component atoms
	for (size_t iq = 0; iq < n_qvals; ++iq) {
		CH_factor->values[iq] =
			table->factors[C].values[iq] + table->factors[H].values[iq];
		CH2_factor->values[iq] =
			table->factors[C].values[iq] + 2 * table->factors[H].values[iq];
		CH3_factor->values[iq] =
			table->factors[C].values[iq] + 3 * table->factors[H].values[iq];
		NH_factor->values[iq] =
			table->factors[N].values[iq] + table->factors[H].values[iq];
		NH2_factor->values[iq] =
			table->factors[N].values[iq] + 2 * table->factors[H].values[iq];
		NH3_factor->values[iq] =
			table->factors[N].values[iq] + 3 * table->factors[H].values[iq];
		OH_factor->values[iq] =
			table->factors[O].values[iq] + table->factors[H].values[iq];
		OH2_factor->values[iq] =
			table->factors[O].values[iq] + 2 * table->factors[H].values[iq];
		SH_factor->values[iq] =
			table->factors[S].values[iq] + table->factors[H].values[iq];
		PO4_factor->values[iq] =
			table->factors[P].values[iq] + 4 * table->factors[O].values[iq];
	}

	// Compute zero, dummy, and vacuum by summing
	CH_factor->zero_ff = table->factors[C].zero_ff + table->factors[H].zero_ff;
	CH2_factor->zero_ff = table->factors[C].zero_ff + 2 * table->factors[H].zero_ff;
	CH3_factor->zero_ff = table->factors[C].zero_ff + 3 * table->factors[H].zero_ff;
	NH_factor->zero_ff = table->factors[N].zero_ff + table->factors[H].zero_ff;
	NH2_factor->zero_ff = table->factors[N].zero_ff + 2 * table->factors[H].zero_ff;
	NH3_factor->zero_ff = table->factors[N].zero_ff + 3 * table->factors[H].zero_ff;
	OH_factor->zero_ff = table->factors[O].zero_ff + table->factors[H].zero_ff;
	OH2_factor->zero_ff = table->factors[O].zero_ff + 2 * table->factors[H].zero_ff;
	SH_factor->zero_ff = table->factors[S].zero_ff + table->factors[H].zero_ff;
	PO4_factor->zero_ff = table->factors[P].zero_ff + 4 * table->factors[O].zero_ff;

	CH_factor->dummy_ff = table->factors[C].dummy_ff + table->factors[H].dummy_ff;
	CH2_factor->dummy_ff = table->factors[C].dummy_ff + 2 * table->factors[H].dummy_ff;
	CH3_factor->dummy_ff = table->factors[C].dummy_ff + 3 * table->factors[H].dummy_ff;
	NH_factor->dummy_ff = table->factors[N].dummy_ff + table->factors[H].dummy_ff;
	NH2_factor->dummy_ff = table->factors[N].dummy_ff + 2 * table->factors[H].dummy_ff;
	NH3_factor->dummy_ff = table->factors[N].dummy_ff + 3 * table->factors[H].dummy_ff;
	OH_factor->dummy_ff = table->factors[O].dummy_ff + table->factors[H].dummy_ff;
	OH2_factor->dummy_ff = table->factors[O].dummy_ff + 2 * table->factors[H].dummy_ff;
	SH_factor->dummy_ff = table->factors[S].dummy_ff + table->factors[H].dummy_ff;
	PO4_factor->dummy_ff = table->factors[P].dummy_ff + 4 * table->factors[O].dummy_ff;

	CH_factor->vacuum_ff =
		table->factors[C].vacuum_ff + table->factors[H].vacuum_ff;
	CH2_factor->vacuum_ff =
		table->factors[C].vacuum_ff + 2 * table->factors[H].vacuum_ff;
	CH3_factor->vacuum_ff =
		table->factors[C].vacuum_ff + 3 * table->factors[H].vacuum_ff;
	NH_factor->vacuum_ff =
		table->factors[N].vacuum_ff + table->factors[H].vacuum_ff;
	NH2_factor->vacuum_ff =
		table->factors[N].vacuum_ff + 2 * table->factors[H].vacuum_ff;
	NH3_factor->vacuum_ff =
		table->factors[N].vacuum_ff + 3 * table->factors[H].vacuum_ff;
	OH_factor->vacuum_ff =
		table->factors[O].vacuum_ff + table->factors[H].vacuum_ff;
	OH2_factor->vacuum_ff =
		table->factors[O].vacuum_ff + 2 * table->factors[H].vacuum_ff;
	SH_factor->vacuum_ff =
		table->factors[S].vacuum_ff + table->factors[H].vacuum_ff;
	PO4_factor->vacuum_ff =
		table->factors[P].vacuum_ff + 4 * table->factors[O].vacuum_ff;
}

void dump_mappings(struct saxs_form_factor_table *table)
{
	struct saxs_resi_atom key = {{0}, {0}};
	enum saxs_ff_type val;
	kh_foreach(table->type_map, key, val,
		   printf("(%s, %s) %d\n",
			  key.residue_name, key.atom_name, val));
}

void print_ff_table(struct saxs_form_factor_table *table)
{
	for (int i = 0; i < HEAVY_ATOM_SIZE; ++i)
	{
		fprintf(stderr, "FFTYPE %2d zero_ff: %f vacuum_ff: %f dummy_ff: %f\n", i,
			table->factors[i].zero_ff,
			table->factors[i].vacuum_ff,
			table->factors[i].dummy_ff);
	}
}

double dummy_ff(const struct saxs_form_factor_table *table,
		const struct mol_atom_group *ag, size_t atom_index)
{
	const struct saxs_form_factor *factor = get_ff(table, ag, atom_index);
	if (factor != NULL) {
		return factor->dummy_ff;
	}

	return 0.0;
}

double vacuum_ff(const struct saxs_form_factor_table *table,
		 const struct mol_atom_group *ag, size_t atom_index)
{
	const struct saxs_form_factor *factor = get_ff(table, ag, atom_index);
	if (factor != NULL) {
		return factor->vacuum_ff;
	}

	return 0.0;
}

void get_dummy_ff(double *dummy_ffs,
                  const struct mol_atom_group *ag, const struct saxs_form_factor_table *table)
{
	for (size_t i = 0; i < ag->natoms; i++) {
		dummy_ffs[i] = dummy_ff(table, ag, i);
	}
}

void get_vacuum_ff(double *vacuum_ffs,
                   const struct mol_atom_group *ag, const struct saxs_form_factor_table *table)
{
	for (size_t i = 0; i < ag->natoms; i++) {
		vacuum_ffs[i] = vacuum_ff(table, ag, i);
	}
}
