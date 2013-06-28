/*
 *	$Id$
 *
 *	Copyright (c) 1991-2013 by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis and F. Wobbe
 *      See LICENSE.TXT file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU Lesser General Public License as published by
 *      the Free Software Foundation; version 3 or any later version.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU Lesser General Public License for more details.
 *
 *      Contact info: www.soest.hawaii.edu/pwessel
 *--------------------------------------------------------------------*/
/* GMT convenience functions used by MATLAB/OCTAVE mex functions
 */
#define STDC_FORMAT_MACROS
#include "gmtmex.h"
#include <math.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <float.h>
#include <math.h>
#include <limits.h>

#if defined(WIN32) && !defined(lrint)
#	define lrint (int64_t)rint
#endif

/* New parser for all GMT mex modules based on design discussed by PW and JL on Mon, 2/21/11 */
/* Wherever we say "Matlab" we mean "Matlab of Octave" */

/* For the Mex interface we will wish to pass either filenames or matrices via GMT command options.
 * We select a Matlab matrix by suppying $ as the file name.  The parser will then find these $
 * arguments and replace them with references to a matrix via the GMT API mechanisms.
 * This requires us to know which options in a module may accept a file name.  As an example,
 * consider surface whose -L option may take a grid.  To pass a Matlab/Octave grid already in memory
 * we would use -L$ and give the grid as an argument to the module, e.g.,
 * Z = surface ('-R0/50/0/50 -I1 -V xyzfile -L$', lowmatrix);
 * For each option that may take a file we need to know what kind of file and if this is input or output.
 * We encode this in a 3-character word, where the first char is the option, the second is the data type,
 * and the third is I(n) or O(out).  E.g., the surface example would have the word LGI.  The data types
 * are P(olygons), L(ines), D(point data), G(rid), C(PT file), T(ext table). [We originally only had
 * D for datasets but perhaps the geometry needs to be passed too (if known); hence the P|L|D char]
 * In addition, the only common option that might take a file is -R which may take a grid as input.
 * We check for that in addition to the module-specific info passed via the key variable.
 *
 * The actual reading/writing will occur in gmt_api.c where we will add a case MEX: for each type.
 * and here we will use mx* for allocation for stuff that is sent to Matlab, and GMT_memory for things
 * that are read and reformatted from Matlab.  This includes packing up GMT grids into Matlab grid structs
 * and vice versa.
 */

#define GMT_MEX_NONE		-3
#define GMT_MEX_EXPLICIT	-2
#define GMT_MEX_IMPLICIT	-1

#define GMT_IS_PS	99	/* Use for PS output; use GMT_IS_GRID or GMT_IS_DATASET for data */


/* Macros for getting the Matlab/Octave ij that correspond to (row,col) [no pad involved] */
/* This one operates on GMT_MATRIX */
#define MEXM_IJ(M,row,col) ((col)*M->n_rows + M->n_rows - (row) - 1)

/* And this on GMT_GRID */
#define MEXG_IJ(M,row,col) ((col)*M->header->ny + M->header->ny - (row) - 1)

int GMTMEX_print_func (FILE *fp, const char *message)
{
	/* Replacement for GMT's gmt_print_func.  It is being used indirectly via
	 * API->print_func.  Purpose of this is to allow Matlab (which cannot use
	 * printf) to reset API->print_func this functions via GMT_Create_Session. */

	mexPrintf (message);
	return 0;
}

int GMTMEX_find_module (void *API, char *module)
{	/* Just search for module and return entry in keys array */
	int k, id = -1;
	for (k = 0; id == -1 && k < N_GMT_MODULES; k++) if (!strcmp (module, module_name[k])) id = k;
	if (id == -1) return (-1);	/* Not found in the known list */
	if ((k = GMT_Call_Module (API, module_name[id], GMT_MODULE_EXIST, NULL)) == GMT_NOERROR) return (id);	/* Found and accessible */
	return (id);	/* Not found in the shared library */
}

char *mxstrdup (const char *s) {
	/* A strdup replacement to be used in Mexs to avoid memory leaks since the Matlab
	   memory management will take care to free the memory allocated by this function */
    char *d = mxMalloc (strlen (s) + 1);
    if (d == NULL) return NULL;
    strcpy (d,s);
    return d;
}

int gmtmex_find_option (char option, char *key[], int n_keys) {
	/* gmtmex_find_option determines if the given option is among the special options that might take $ as filename */
	int pos = -1, k;
	for (k = 0; pos == -1 && k < n_keys; k++) if (key[k][0] == option) pos = k;
	return (pos);	/* -1 if not found, otherwise the position in the key array */
}

int gmtmex_get_arg_pos (char *arg)
{	/* Look for a $ in the arg; if found return position, else return -1. Skips $ inside quoted texts */
	int pos, k;
	unsigned int mute = 0;
	for (k = 0, pos = -1; pos == -1 && k < (int)strlen (arg); k++) {
		if (arg[k] == '\"' || arg[k] == '\'') mute = !mute;	/* Do not consider $ inside quotes */
		if (!mute && arg[k] == '$') pos = k;	/* Found a $ sign */
	}
	return (pos);	/* Either -1 (not found) or in the 0-(strlen(arg)-1) range [position of $] */
}

unsigned int gmtmex_get_key_pos (char *key[], unsigned int n_keys, struct GMT_OPTION *head, int def[])
{	/* Must determine if default input and output have been set via program options or if they should be added explicitly.
 	 * As an example, consider the GMT command grdfilter in.nc -Fg200k -Gfilt.nc.  In Matlab this might be
	 * filt = GMT_grdfilter ('$ -Fg200k -G$', in);
	 * However, it is more natural not to specify the lame -G$, i.e.
	 * filt = GMT_grdfilter ('$ -Fg200k', in);
	 * In that case we need to know that -G is the default way to specify the output grid and if -G is not given we
	 * must associate -G with the first left-hand-side item (here filt).
	 */
	int pos, PS = 0;
	struct GMT_OPTION *opt = NULL;
	def[GMT_IN] = def[GMT_OUT] = GMT_MEX_IMPLICIT;	/* Initialize to setting the i/o implicitly */
	
	/* Loop over the module options to see if inputs and outputs are set explicitly or implicitly */
	for (opt = head; opt; opt = opt->next) {
		pos = gmtmex_find_option (opt->option, key, n_keys);	/* First see if this option is one that might take $ */
		if (pos == -1) continue;		/* No, it was some other harmless option, e.g., -J, -O ,etc. */
		/* Here, the current option is one that might take an input or output file. See if it matches
		 * the UPPERCASE I or O [default source/dest] rather than the standard i|o (optional input/output) */
		if (key[pos][2] == 'I')  /* Default input  is actually set explicitly via option setting now indicated by key[pos] */
			def[GMT_IN] = GMT_MEX_EXPLICIT;
		else if (key[pos][2] == 'O')  /* Default output is actually set explicitly via option setting now indicated by key[pos] */
			def[GMT_OUT] = GMT_MEX_EXPLICIT;
	}
	/* Here, if def[] == GMT_MEX_IMPLICIT (the default in/out option was NOT given), then we want to return the corresponding entry in key */
	for (pos = 0; pos < n_keys; pos++) {	/* For all module options that might take a file */
		if ((key[pos][2] == 'I' || key[pos][2] == 'i') && key[pos][0] == '-') /* This program takes no input (e.g., psbasemap) */
			def[GMT_IN]  = GMT_MEX_NONE;
		else if (key[pos][2] == 'I' && def[GMT_IN]  == GMT_MEX_IMPLICIT) /* Must add implicit input; use def to determine option,type */
			def[GMT_IN]  = pos;
		if ((key[pos][2] == 'O' || key[pos][2] == 'o') && key[pos][0] == '-') /* This program produces no output */
			def[GMT_OUT] = GMT_MEX_NONE;
		else if (key[pos][2] == 'O' && def[GMT_OUT] == GMT_MEX_IMPLICIT) /* Must add implicit output; use def to determine option,type */
			def[GMT_OUT] = pos;
		if ((key[pos][2] == 'O' || key[pos][2] == 'o') && key[pos][1] == 'X' && key[pos][0] == '-') PS = 1;	/* This program produces PostScript */
	}
	return (PS);
}

int gmtmex_get_arg_dir (char option, char *key[], int n_keys, int *data_type, int *geometry)
{
	int item;
	
	/* 1. First determine if this option is one of the choices in key */
	
	item = gmtmex_find_option (option, key, n_keys);
	if (item == -1)		/* This means a coding error we must fix */
		mexErrMsgTxt ("GMTMEX_pre_process: This option does not allow $ arguments\n");
	
	/* 2. Assign direction, data_type, and geometry */
	
	switch (key[item][1]) {	/* 2nd char contains the data type code */
		case 'G':
			*data_type = GMT_IS_GRID;
			*geometry = GMT_IS_SURFACE;
			break;
		case 'P':
			*data_type = GMT_IS_DATASET;
			*geometry = GMT_IS_POLY;
			break;
		case 'L':
			*data_type = GMT_IS_DATASET;
			*geometry = GMT_IS_LINE;
			break;
		case 'D':
			*data_type = GMT_IS_DATASET;
			*geometry = GMT_IS_POINT;
			break;
		case 'C':
			*data_type = GMT_IS_CPT;
			*geometry = GMT_IS_NONE;
			break;
		case 'T':
			*data_type = GMT_IS_TEXTSET;
			*geometry = GMT_IS_NONE;
			break;
		case 'I':
			*data_type = GMT_IS_IMAGE;
			*geometry = GMT_IS_SURFACE;
			break;
		case 'X':
			*data_type = GMT_IS_PS;
			*geometry = GMT_IS_NONE;
			break;
		default:
			mexErrMsgTxt ("GMTMEX_pre_process: Bad data_type character in 3-char module code!\n");
			break;
	}
	/* Third key character contains the in/out code */
	if (key[item][2] == 'I') key[item][2] = 'i';	/* This was the default input option set explicitly; no need to add later */
	if (key[item][2] == 'O') key[item][2] = 'o';	/* This was the default output option set explicitly; no need to add later */
	return ((key[item][2] == 'i') ? GMT_IN : GMT_OUT);
}

char **make_char_array (char *string, unsigned int *n_items, char type)
{	/* Turn the comma-separated list of 3-char codes into an array of such codes.
 	 * In the process, replace any ?-types with type. */
	size_t len, k, n;
	char **s = NULL;
	char *next, *tmp;
	
	if (!string) return NULL;
	len = strlen (string);
	if (len == 0) return NULL;
	tmp = strdup (string);
	if (type) for (k = 0; k < strlen (tmp); k++) if (tmp[k] == '?') tmp[k] = type;
	for (k = n = 0; k < len; k++) if (tmp[k] == ',') n++;
	n++;
	s = (char **) calloc (n, sizeof (char *));
	k = 0;
	while ((next = strsep (&tmp, ",")) != NULL) {
		s[k++] = strdup (next);
	}
	*n_items = (unsigned int)n;
	free ((void *)tmp);
	return s;
}

struct GMT_GRID *GMTMEX_grid_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Grid container and associate it with a Matlab grid.
 	 * If direction is GMT_IN then we are given a Matlab grid and can determine size etc. */
	struct GMT_GRID *G = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		mxArray *mx_ptr = NULL;
		double *inc = NULL, *range = NULL, *reg = NULL;
		unsigned int registration;
		
		if (mxIsEmpty (ptr)) mexErrMsgTxt ("The input that was supposed to contain the Grid, is empty\n");
		if (!mxIsStruct (ptr)) mexErrMsgTxt ("Expected a Grid for input\n");
		mx_ptr = mxGetField (ptr, 0, "inc");
		if (mx_ptr == NULL) mexErrMsgTxt ("Could not find inc array with Grid increments\n");
		inc = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "range");
		if (mx_ptr == NULL) mexErrMsgTxt ("Could not find range array for Grid range\n");
		range = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "registration");
		if (mx_ptr == NULL) mexErrMsgTxt ("Could not find registration array for Grid registration\n");
		reg = mxGetData (mx_ptr);
		registration = (unsigned int)lrint (reg[0]);
		if ((G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, NULL, range, inc, registration, 0, NULL)) == NULL)
			mexErrMsgTxt ("Failure to alloc GMT source matrix\n");
		mx_ptr = mxGetField (ptr, 0, "data");
		if (mx_ptr == NULL) mexErrMsgTxt ("Could not find data array for Grid\n");
		G->data = mxGetData (mx_ptr);
		G->alloc_mode = GMT_ALLOCATED_EXTERNALLY;	/* Since array was allocated by Matlab */
		GMT_Report (API, GMT_MSG_DEBUG, " Allocate GMT Grid %lx in gmtmex_parser\n", (long)G);
		GMT_Report (API, GMT_MSG_DEBUG, " Registered GMT Grid array %lx via memory reference from Matlab\n", (long)G->data);
	}
	else {	/* Set dummy range and inc so that GMT_check_region in GMT_Register_IO will pass */
		double range[4] = {0.0, 0.0, 0.0, 0.0}, inc[2] = {1.0, 1.0};
		unsigned int registration = GMT_GRID_NODE_REG;
		if ((G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY | GMT_VIA_OUTPUT, NULL, range, inc, registration, 0, NULL)) == NULL)
			mexErrMsgTxt ("Failure to alloc GMT source matrix\n");
	}
	return (G);
}

struct GMT_MATRIX *GMTMEX_matrix_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Matrix container and associate it with a Matlab matrix.
 	 * If direction is GMT_IN then we are given a Matlab matrix and can determine size etc.
	 * If output then we dont know size but we can specify type */
	uint64_t dim[2] = {0, 0};
	unsigned int mode = 0;
	struct GMT_MATRIX *M = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known */
		if (!mxIsNumeric (ptr)) mexErrMsgTxt ("Expected a Matrix for input\n");
		dim[0] = mxGetN (ptr);
		dim[1] = mxGetM (ptr);
	}
	else
		mode = GMT_VIA_OUTPUT;
	if ((M = GMT_Create_Data (API, GMT_IS_MATRIX, GMT_IS_SURFACE, mode, dim, NULL, NULL, 0, 0, NULL)) == NULL)
		mexErrMsgTxt ("Failure to alloc GMT source matrix\n");

	GMT_Report (API, GMT_MSG_DEBUG, " Allocate GMT Matrix %lx in gmtmex_parser\n", (long)M);
	M->n_rows    = dim[1];
	M->n_columns = dim[0];
	/* Set a silly range so that GMT_check_region in GMT_Register_IO will pass */
	M->range[GMT_XLO] = M->range[GMT_YLO] = 1.0;
	M->range[GMT_XHI] = (double)M->n_columns;	M->range[GMT_YHI] = (double)M->n_rows;
	if (direction == GMT_IN) {	/* We can inquire about the input */
		if (mxIsDouble(ptr)) {
			M->type = GMT_DOUBLE;
			M->data.f8 = mxGetData (ptr);
		}
		else if (mxIsSingle(ptr)) {
			M->type = GMT_FLOAT;
			M->data.f4 = (float *)mxGetData (ptr);
		}
		else if (mxIsInt32(ptr)) {
			M->type = GMT_INT;
			M->data.si4 = (int32_t *)mxGetData (ptr);
		}
		else if (mxIsInt16(ptr)) {
			M->type = GMT_SHORT;
			M->data.si2 = (int16_t *)mxGetData (ptr);
		}
		else if (mxIsInt8(ptr)) {
			M->type = GMT_CHAR;
			M->data.sc1 = (int8_t *)mxGetData (ptr);
		}
		else
			mexErrMsgTxt ("Unsupported data type in GMT matrix input.");

		M->shape = GMT_IS_COL_FORMAT;
		M->dim = M->n_rows;	// This is actualy wrong if input data is scanline as for Octave oct
		M->alloc_mode = GMT_ALLOCATED_EXTERNALLY;	/* Since matrix was allocated by Matlab */
	}
	else {	/* On output we produce doubles */
		M->type = GMT_FLOAT;
	}
	return (M);
}

int GMTMEX_Register_IO (void *API, unsigned int data_type, unsigned int geometry, unsigned int direction, const mxArray *ptr)
{	/* Create the grid or matrix contains, register them, and return the ID */
	int ID = GMT_NOTSET;
	struct GMT_GRID *G = NULL;		/* Pointer to grid container */
	struct GMT_MATRIX *M = NULL;		/* Pointer to matrix container */

	switch (data_type) {
		case GMT_IS_GRID:
			G = GMTMEX_grid_init (API, direction, ptr);	/* Get a grid and associate it with the Matlab grid pointer (if input) */
			ID = GMT_Get_ID (API, GMT_IS_GRID, direction, G);
			break;
		case GMT_IS_DATASET:
			M = GMTMEX_matrix_init (API, direction, ptr);	/* Get a matrix container and associate it with the Matlab pointer (if input) */
			if ((ID = GMT_Register_IO (API, GMT_IS_DATASET, GMT_IS_REFERENCE + GMT_VIA_MATRIX, geometry, direction, NULL, M)) == GMT_NOTSET) {
				mexErrMsgTxt ("GMTMEX_pre_process: Failure to register GMT matrix source or destination\n");
			}
			break;
		default:
			mexErrMsgTxt ("GMTMEX_pre_process: Bad data type\n");
			break;
	}
	return (ID);
}

int GMTMEX_pre_process (void *API, const char *module, mxArray *plhs[], int nlhs, const mxArray *prhs[], int nrhs, char *keys, struct GMT_OPTION *head, struct GMTMEX **X)
{
	/* API controls all things within GMT.
	 * plhs (and nlhs) are the outputs specified on the left side of the equal sign in Matlab.
	 * prhs (and nrhs) are the inputs specified after the option string in the GMT-mex function.
	 * keys is comma-separated string with 3-char codes for current module i/o.
	 * opt is the linked list of GMT options passed in.
	 */
	
	int lr_pos[2] = {0, 0};	/* These position keeps track where we are in the L and R pointer arrays */
	int direction;		/* Either GMT_IN or GMT_OUT */
	int data_type;		/* Either GMT_IS_DATASET, GMT_IS_TEXTSET, GMT_IS_GRID, GMT_IS_CPT, GMT_IS_IMAGE */
	int geometry;		/* Either GMT_IS_NONE, GMT_IS_POINT, GMT_IS_LINE, GMT_IS_POLY, or GMT_IS_SURFACE */
	int def[2];		/* Either GMT_MEX_EXPLICIT or the item number in the keys array */
	int ID, error;
	int pos;
	unsigned int k, n_keys = 0, PS, n_alloc = 8U, n_items = 0;
	char name[GMT_STR16];	/* Used to hold the GMT API embedded file name, e.g., @GMTAPI@-###### */
	char **key = NULL;
	char *text = NULL;
	char type = 0;
	void *ptr = NULL;	
	struct GMT_OPTION *opt, *new_ptr;	/* Pointer to a GMT option structure */
	struct GMT_GRID *G = NULL;		/* Pointer to grid container */
	struct GMT_MATRIX *M = NULL;		/* Pointer to matrix container */
	struct GMTMEX *info = NULL;

	if (!strcmp (module, "gmtread") || !strcmp (module, "gmtwrite"))  {	/* Special case: Must determine which data type we are dealing with */
		struct GMT_OPTION *t_ptr;
		if ((t_ptr = GMT_Find_Option (API, 'T', head))) {	/* Found the -T<type> option */
			type = toupper (t_ptr->arg[0]);	/* Find type and replace ? in keys with this type in uppercase (DGCIT) in make_char_array below */
		}
		if (!strcmp (module, "gmtwrite") && (t_ptr = GMT_Find_Option (API, GMT_OPT_INFILE, head))) {	/* Found a -<<file> option; this is actually the output file */
			t_ptr->option = GMT_OPT_OUTFILE;
		}
	}
	
	key = make_char_array (keys, &n_keys, type);
	info = malloc (n_alloc * sizeof (struct GMTMEX));
	
	PS = gmtmex_get_key_pos (key, n_keys, head, def);	/* Determine if we must add the primary in and out arguments to the option list */
	for (direction = GMT_IN; direction <= GMT_OUT; direction++) {
		if (def[direction] == GMT_MEX_NONE) continue;	/* No source or destination required */
		if (def[direction] == GMT_MEX_EXPLICIT) continue;	/* Source or destination was set explicitly; skip */
		/* Must add the primary input or output from prhs[0] or plhs[0] */
		(void)gmtmex_get_arg_dir (key[def[direction]][0], key, n_keys, &data_type, &geometry);		/* Get info about the data set */
		ptr = (direction == GMT_IN) ? (void *)prhs[lr_pos[direction]] : (void *)plhs[lr_pos[direction]];	/* Pick the next left or right side pointer */
		ID = GMTMEX_Register_IO (API, data_type, geometry, direction, ptr);
		/* Register a Matlab/Octave entity as a source or destination */
		if (direction == GMT_OUT) {
			if (n_items == n_alloc) info = realloc ((void *)info, (n_alloc += 8) * sizeof (struct GMTMEX));
			info[n_items].type = data_type;
			info[n_items].ID = ID;
			info[n_items].lhs_index = lr_pos[GMT_OUT];
			n_items++;
		}
		lr_pos[direction]++;		/* Advance position counter for next time */
		if (GMT_Encode_ID (API, name, ID) != GMT_NOERROR) {	/* Make filename with embedded object ID */
			mexErrMsgTxt ("GMTMEX_pre_process: Failure to encode string\n");
		}
		new_ptr = GMT_Make_Option (API, key[def[direction]][0], name);	/* Create the missing (implicit) GMT option */
		GMT_Append_Option (API, new_ptr, head);				/* Append it to the option list */
	}
		
	for (opt = head; opt; opt = opt->next) {	/* Loop over the module options given */
		if (PS && opt->option == GMT_OPT_OUTFILE) PS++;
		/* Determine if this option has a $ in its argument and if so return its position in pos; return -1 otherwise */
		if ((pos = gmtmex_get_arg_pos (opt->arg)) == -1) continue;	/* No $ argument found or it is part of a text string */
		
		/* Determine several things about this option, such as direction, data type, and geometry */
		direction = gmtmex_get_arg_dir (opt->option, key, n_keys, &data_type, &geometry);
		ptr = (direction == GMT_IN) ? (void *)prhs[lr_pos[direction]] : (void *)plhs[lr_pos[direction]];	/* Pick the next left or right side pointer */
		ID = GMTMEX_Register_IO (API, data_type, geometry, direction, ptr);
		if (direction == GMT_OUT) {
			if (n_items == n_alloc) info = realloc ((void *)info, (n_alloc += 8) * sizeof (struct GMTMEX));
			info[n_items].type = data_type;
			info[n_items].ID = ID;
			info[n_items].lhs_index = lr_pos[GMT_OUT];
			n_items++;
		}
		if (GMT_Encode_ID (API, name, ID) != GMT_NOERROR) {	/* Make filename with embedded object ID */
			mexErrMsgTxt ("GMTMEX_pre_process: Failure to encode string\n");
		}
		lr_pos[direction]++;		/* Advance position counter for next time */
		
		/* Replace the option argument with the embedded file */
		if (GMT_Update_Option (API, opt, name)) {
			mexErrMsgTxt ("GMTMEX_pre_process: Failure to update option argument\n");
		}
	}
	if (n_items && n_items < n_alloc) info = realloc ((void *)info, n_items * sizeof (struct GMTMEX));
	else if (n_items == 0) free ((void *)info);
	
	if (PS == 1)	/* No redirection of PS to a file */
		error = 1;
	else if (PS > 2)	/* Too many output files for PS */
		error = 2;
	else
		error = GMT_NOERROR;
	for (k = 0; k < n_keys; k++) free ((void *)key[k]);
	free ((void *)key);
	
	text = GMT_Create_Cmd (API, head);
	GMT_Report (API, GMT_MSG_VERBOSE, "Args are now [%s]\n", text);
	GMT_Destroy_Cmd (API, &text);

	/* Here, a command line '-F200k -G$ $ -L$ -P' has been changed to '-F200k -G@GMTAPI@-000001 @GMTAPI@-000002 -L@GMTAPI@-000003 -P'
	 * where the @GMTAPI@-00000x are encodings to registered resources or destinations */
	*X = info;
	return (error ? -error : n_items);
}

int GMTMEX_post_process (void *API, struct GMTMEX *X, int n_items, mxArray *plhs[])
{	/* Get the data from GMT output items into the corresponding Matlab struct or matrix */
	int item, k, n;
	unsigned int row, col;
	uint64_t gmt_ij;
	float  *f = NULL;
	double *d = NULL, *dptr = NULL, *G_x = NULL, *G_y = NULL, *x = NULL, *y = NULL;
	struct GMT_GRID *G = NULL;
	struct GMT_MATRIX *M = NULL;
	mxArray *mxGrd = NULL, *mx_x = NULL, *mx_y= NULL;
	mxArray *mxProjectionRef = NULL;
	mxArray *mxHeader = NULL, *mxtmp = NULL;
	mxArray *grid_struct = NULL;
	char    *fieldnames[21];	/* this array contains the names of the fields of the output structure. */
	
	GMT_Report (API, GMT_MSG_VERBOSE, "Enter GMTMEX_post_process\n");
	for (item = 0; item < n_items; item++) {
		k = X[item].lhs_index;	/* Short-hand for index into plhs[] */
		switch (X[item].type) {
			case GMT_IS_GRID:	/* Return grids via a float (mxSINGLE_CLASS) matrix in a struct */
				if ((G = GMT_Retrieve_Data (API, X[item].ID)) == NULL) mexErrMsgTxt ("Error retrieving grid from GMT\n");
				/* Create a Matlab struct for this grid */
				fieldnames[0]  = mxstrdup ("ProjectionRefPROJ4");
				fieldnames[1]  = mxstrdup ("ProjectionRefWKT");
				fieldnames[2]  = mxstrdup ("hdr");
				fieldnames[3]  = mxstrdup ("range");
				fieldnames[4]  = mxstrdup ("inc");
				fieldnames[5]  = mxstrdup ("dim");
				fieldnames[6]  = mxstrdup ("n_rows");
				fieldnames[7]  = mxstrdup ("n_columns");
				fieldnames[8]  = mxstrdup ("MinMax");
				fieldnames[9]  = mxstrdup ("NoDataValue");
				fieldnames[10] = mxstrdup ("registration");
				fieldnames[11] = mxstrdup ("title");
				fieldnames[12] = mxstrdup ("remark");
				fieldnames[13] = mxstrdup ("command");
				fieldnames[14] = mxstrdup ("DataType");
				fieldnames[15] = mxstrdup ("LayerCount");
				fieldnames[16] = mxstrdup ("ColorMap");
				fieldnames[17] = mxstrdup ("data");
				fieldnames[18] = mxstrdup ("x_units");
				fieldnames[19] = mxstrdup ("y_units");
				fieldnames[20] = mxstrdup ("z_units");
				grid_struct = mxCreateStructMatrix (1, 1, 21, (const char **)fieldnames );
				mxHeader    = mxCreateNumericMatrix (1, 9, mxDOUBLE_CLASS, mxREAL);
				dptr = mxGetPr (mxHeader);
				dptr[0] = 0;	dptr[1] = 1;	dptr[2] = 0;	dptr[3] = 2;
				mxSetField (grid_struct, 0, "hdr", mxHeader);

				dptr = mxGetPr(mxtmp = mxCreateNumericMatrix (1, 4, mxDOUBLE_CLASS, mxREAL));
				for (n = 0; n < 4; n++) dptr[n] = G->header->wesn[n];
				mxSetField (grid_struct, 0, "range", mxtmp);

				dptr = mxGetPr(mxtmp = mxCreateNumericMatrix (1, 2, mxDOUBLE_CLASS, mxREAL));
				for (n = 0; n < 2; n++) dptr[n] = G->header->inc[n];
				mxSetField (grid_struct, 0, "inc", mxtmp);

				mxtmp = mxCreateDoubleScalar ((double)G->header->ny); 
				mxSetField (grid_struct, 0, (const char *) "n_rows", mxtmp);

				mxtmp = mxCreateDoubleScalar ((double)G->header->nx); 
				mxSetField (grid_struct, 0, (const char *) "n_columns", mxtmp);

				dptr = mxGetPr(mxtmp = mxCreateNumericMatrix (1, 2, mxDOUBLE_CLASS, mxREAL));
				dptr[0] = G->header->z_min;	dptr[1] = G->header->z_max;
				mxSetField (grid_struct, 0, "MinMax", mxtmp);

				mxtmp = mxCreateDoubleScalar ((double)G->header->nan_value); 
				mxSetField (grid_struct, 0, (const char *) "NoDataValue", mxtmp);

				dptr = mxGetPr(mxtmp = mxCreateNumericMatrix (1, 2, mxDOUBLE_CLASS, mxREAL));
				dptr[0] = G->header->ny;	dptr[1] = G->header->nx;
				mxSetField (grid_struct, 0, "dim", mxtmp);

				mxtmp = mxCreateDoubleScalar ((double)G->header->registration); 
				mxSetField (grid_struct, 0, (const char *) "registration", mxtmp);

				mxtmp = mxCreateString (G->header->title);
				mxSetField (grid_struct, 0, (const char *) "title", mxtmp);

				mxtmp = mxCreateString (G->header->command);
				mxSetField (grid_struct, 0, (const char *) "command", mxtmp);

				mxtmp = mxCreateString (G->header->remark);
				mxSetField (grid_struct, 0, (const char *) "remark", mxtmp);

				mxtmp = mxCreateString (G->header->x_units);
				mxSetField (grid_struct, 0, (const char *) "x_units", mxtmp);

				mxtmp = mxCreateString (G->header->y_units);
				mxSetField (grid_struct, 0, (const char *) "y_units", mxtmp);

				mxtmp = mxCreateString (G->header->z_units);
				mxSetField (grid_struct, 0, (const char *) "z_units", mxtmp);

				if (!G->data)
					mexErrMsgTxt ("GMTMEX_post_process: programming error, output matrix G is empty\n");

				mxGrd = mxCreateNumericMatrix (G->header->ny, G->header->nx, mxSINGLE_CLASS, mxREAL);
				f = mxGetData (mxGrd);
				/* Load the real grd array into a double matlab array by transposing 
			           from unpadded GMT grd format to unpadded matlab format */
				for (gmt_ij = row = 0; row < G->header->ny; row++) for (col = 0; col < G->header->nx; col++, gmt_ij++)
					f[MEXG_IJ(G,row,col)] = G->data[gmt_ij];
				mxSetField (grid_struct, 0, "z", mxGrd);
				/* Also return x and y arrays */
				G_x = GMT_Get_Coord (API, GMT_IS_GRID, GMT_X, G);	/* Get array of x coordinates */
				G_y = GMT_Get_Coord (API, GMT_IS_GRID, GMT_Y, G);	/* Get array of y coordinates */
				mx_x = mxCreateNumericMatrix (1, G->header->nx, mxDOUBLE_CLASS, mxREAL);
				mx_y = mxCreateNumericMatrix (1, G->header->ny, mxDOUBLE_CLASS, mxREAL);
				x = mxGetData (mx_x);
				y = mxGetData (mx_y);
				memcpy (x, G_x, G->header->nx * sizeof (double));
				memcpy (y, G_y, G->header->ny * sizeof (double));
				GMT_Destroy_Data (API, G_x);
				GMT_Destroy_Data (API, G_y);
				mxSetField (grid_struct, 0, "x", mx_x);
				mxSetField (grid_struct, 0, "y", mx_y);
				plhs[k] = grid_struct;
				break;
				
			case GMT_IS_DATASET:	/* Return tables with double (mxDOUBLE_CLASS) matrix */
				if ((M = GMT_Retrieve_Data (API, X[item].ID)) == NULL) mexErrMsgTxt ("Error retrieving matrix from GMT\n");
				/* Create a Matlab matrix to hold this GMT matrix */
				plhs[k] = mxCreateNumericMatrix (M->n_rows, M->n_columns, mxDOUBLE_CLASS, mxREAL);
				d = mxGetData (plhs[k]);
				/* Load the real data matrix into a double matlab array copying columns */
				for (col = 0; col < M->n_columns; col++)
					memcpy (&d[M->n_rows*col], M->data.f8, M->n_rows * sizeof(double));
				break;
		}
	}
	GMT_Report (API, GMT_MSG_VERBOSE, "Exit GMTMEX_post_process\n");
	return (0);	/* Maybe we should turn this function to void but the gmt5 wrapper expects a return value */
}
