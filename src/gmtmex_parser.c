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

#include "gmt_mex.h"

//#define TESTING
#ifdef TESTING
unsigned int unique_ID = 0;
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

int GMTMEX_print_func (FILE *fp, const char *message)
{
	/* Replacement for GMT's gmt_print_func.  It is being used indirectly via
	 * API->print_func.  Purpose of this is to allow Matlab (which cannot use
	 * printf) to reset API->print_func this functions via GMT_Create_Session. */

	mexPrintf (message);
	return 0;
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
	for (k = 0, pos = -1; pos == -1 && k < strlen (arg); k++) {
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
	
	for (opt = head; opt; opt = opt->next) {	/* Loop over the module options to see if inputs and outputs are set explicitly or implicitly */
		pos = gmtmex_find_option (opt->option, key, n_keys);	/* First see if this option is one that might take $ */
		if (pos == -1) continue;		/* No, it was some other harmless option, e.g., -J, -O ,etc. */
		/* Here, the current option is one that might take an input or output file. See if it matches
		 * the UPPERCASE I or O [default source/dest] rather than the standard i|o (optional input/output) */
		if (key[pos][2] == 'I') def[GMT_IN]  = GMT_MEX_EXPLICIT;	/* Default input  is actually set explicitly via option setting now indicated by key[pos] */
		if (key[pos][2] == 'O') def[GMT_OUT] = GMT_MEX_EXPLICIT;	/* Default output is actually set explicitly via option setting now indicated by key[pos] */
	}
	/* Here, if def[] == GMT_MEX_IMPLICIT (the default in/out option was NOT given), then we want to return the corresponding entry in key */
	for (pos = 0; pos < n_keys; pos++) {	/* For all module options that might take a file */
		if ((key[pos][2] == 'I' || key[pos][2] == 'i') && key[pos][0] == '-') def[GMT_IN]  = GMT_MEX_NONE;	/* This program takes no input (e.g., psbasemap) */
		else if (key[pos][2] == 'I' && def[GMT_IN]  == GMT_MEX_IMPLICIT) def[GMT_IN]  = pos;	/* Must add implicit input; use def to determine option,type */
		if ((key[pos][2] == 'O' || key[pos][2] == 'o') && key[pos][0] == '-') def[GMT_OUT] = GMT_MEX_NONE;	/* This program produces no output */
		else if (key[pos][2] == 'O' && def[GMT_OUT] == GMT_MEX_IMPLICIT) def[GMT_OUT] = pos;	/* Must add implicit output; use def to determine option,type */
		if ((key[pos][2] == 'O' || key[pos][2] == 'o') && key[pos][1] == 'X' && key[pos][0] == '-') PS = 1;	/* This program produces PostScript */
	}
	return (PS);
}

int gmtmex_get_arg_dir (char option, char *key[], int n_keys, int *data_type, int *geometry)
{
	int item;
	
	/* 1. First determine if this option is one of the choices in key */
	
	item = gmtmex_find_option (option, key, n_keys);
	if (item == -1) fprintf (stderr, "GMTMEX_pre_process: This option does not allow $ arguments\n");	/* This means a coding error we must fix */
	
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
			fprintf (stderr, "GMTMEX_pre_process: Bad data_type character in 3-char module code!\n");
			break;
	}
	/* Third key character contains the in/out code */
	if (key[item][2] == 'I') key[item][2] = 'i';	/* This was the default input option set explicitly; no need to add later */
	if (key[item][2] == 'O') key[item][2] = 'o';	/* This was the default output option set explicitly; no need to add later */
	return ((key[item][2] == 'i') ? GMT_IN : GMT_OUT);
}

char ** make_char_array (char *string, unsigned int *n_items)
{
	unsigned int len, k, n;
	char **s = NULL;
	char *next, *tmp;
	
	if (!string) return NULL;
	len = strlen (string);
	if (len == 0) return NULL;
	tmp = strdup (string);
	for (k = n = 0; k < len; k++) if (tmp[k] == ',') n++;
	n++;
	s = (char **) calloc (n, sizeof (char *));
	k = 0;
	while ((next = strsep (&tmp, ",")) != NULL) {
		s[k++] = strdup (next);
	}
	*n_items = n;
	free ((void *)tmp);
	return s;
}

struct GMT_MATRIX *GMTMEX_matrix_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Matrix container and associate it with a Matlab matrix.
 	 * If direction is GMT_IN then we are given a Matlab matrix and can determine size etc.
	 * If output then we dont know size but we can specify type */
	int col, in_ID;
	uint64_t dim[2] = {0, 0};
	struct GMT_MATRIX *M = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known */
		dim[0] = mxGetN (ptr);
		dim[1] = mxGetM (ptr);
	}
	if ((M = GMT_Create_Data (API, GMT_IS_MATRIX, GMT_IS_SURFACE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
		mexErrMsgTxt ("Failure to alloc GMT source matrix\n");

	M->n_rows = dim[1];
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
	}
	else {	/* On output we produce doubles */
		M->type = GMT_FLOAT;
		M->alloc_mode = GMT_REFERENCE;
	}
	return (M);
}

int GMTMEX_pre_process (void *API, mxArray *plhs[], int nlhs, const mxArray *prhs[], int nrhs, char *keys, struct GMT_OPTION *head, struct GMTMEX **X)
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
	unsigned int k, n_keys = 0, pos, PS, n_alloc = 8U, n_items = 0;
	char name[GMT_STR16];	/* Used to hold the GMT API embedded file name, e.g., @GMTAPI@-###### */
	char **key = NULL;
	struct GMT_OPTION *opt, *new_ptr;	/* Pointer to a GMT option structure */
	struct GMT_MATRIX *M = NULL;		/* Pointer to matrix container */
	struct GMTMEX *info = NULL;
#ifndef TESTING
	void *ptr = NULL;		/* Void pointer used to point to either L or R side pointer argument */
#endif
	
	key = make_char_array (keys, &n_keys);
	info = malloc (n_alloc * sizeof (struct GMTMEX));
	
	PS = gmtmex_get_key_pos (key, n_keys, head, def);	/* Determine if we must add the primary in and out arguments to the option list */
	for (direction = GMT_IN; direction <= GMT_OUT; direction++) {
		if (def[direction] == GMT_MEX_NONE) continue;	/* No source or destination required */
		if (def[direction] == GMT_MEX_EXPLICIT) continue;	/* Source or destination was set explicitly; skip */
		/* Must add the primary input or output from prhs[0] or plhs[0] */
		(void)gmtmex_get_arg_dir (key[def[direction]][0], key, n_keys, &data_type, &geometry);		/* Get info about the data set */
#ifdef TESTING
		ID = unique_ID++;
#else
		ptr = (direction == GMT_IN) ? (void *)prhs[lr_pos[direction]] : (void *)plhs[lr_pos[direction]];	/* Pick the next left or right side pointer */
		M = GMTMEX_matrix_init (API, direction, ptr);	/* Get a matrix container and associate it with the Matlab pointer */
		/* Register a Matlab/Octave entity as a source or destination */
		if ((ID = GMT_Register_IO (API, data_type, GMT_IS_REFERENCE + GMT_VIA_MATRIX, geometry, direction, NULL, M)) == GMT_NOTSET) {
			mexErrMsgTxt ("GMTMEX_pre_process: Failure to register GMT source or destination\n");
		}
		M = NULL;	/* Just to be nice and clean */
		if (direction == GMT_OUT) {
			if (n_items == n_alloc) info = realloc ((void *)info, (n_alloc += 8) * sizeof (struct GMTMEX));
			info[n_items].type = data_type;
			info[n_items].ID = ID;
			info[n_items].lhs_index = lr_pos[direction];
			n_items++;
		}
#endif
		lr_pos[direction]++;		/* Advance position counter for next time */
		if (GMT_Encode_ID (API, name, ID) != GMT_NOERROR) {	/* Make filename with embedded object ID */
			mexErrMsgTxt ("GMTMEX_pre_process: Failure to encode string\n");
		}
		new_ptr = GMT_Make_Option (API, key[def[direction]][0], name);	/* Create the missing (implicit) GMT option */
		GMT_Append_Option (API, new_ptr, head);				/* Append it to the option list */
	}
		
	for (opt = head; opt; opt = opt->next) {	/* Loop over the module options given */
		if (PS && opt->option == GMT_OPT_OUTFILE) PS++;
		/* Determine if this option as a $ in its argument and if so return its position in pos; return -1 otherwise */
		if ((pos = gmtmex_get_arg_pos (opt->arg)) == -1) continue;	/* No $ argument found or it is part of a text string */
		
		/* Determine several things about this option, such as direction, data type, method, and geometry */
		direction = gmtmex_get_arg_dir (opt->option, key, n_keys, &data_type, &geometry);
#ifdef TESTING
		ID = unique_ID++;
#else
		ptr = (direction == GMT_IN) ? (void *)prhs[lr_pos[direction]] : (void *)plhs[lr_pos[direction]];	/* Pick the next left or right side pointer */
		M = GMTMEX_matrix_init (API, direction, ptr);	/* Get a matrix container and associate it with the Matlab pointer */
		/* Register a Matlab/Octave entity as a source or destination */
		if ((ID = GMT_Register_IO (API, data_type, GMT_IS_REFERENCE + GMT_VIA_MATRIX, geometry, direction, NULL, M)) == GMT_NOTSET) {
			mexErrMsgTxt ("GMTMEX_pre_process: Failure to register GMT source or destination\n");
		}
		M = NULL;	/* Just to be nice and clean */
		if (direction == GMT_OUT) {
			if (n_items == n_alloc) info = realloc ((void *)info, (n_alloc += 8) * sizeof (struct GMTMEX));
			info[n_items].type = data_type;
			info[n_items].ID = ID;
			info[n_items].lhs_index = lr_pos[direction];
			n_items++;
		}
#endif
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
	
	/* Here, a command line '-F200k -G$ $ -L$ -P' has been changed to '-F200k -G@GMTAPI@-000001 @GMTAPI@-000002 -L@GMTAPI@-000003 -P'
	 * where the @GMTAPI@-00000x are encodings to registered resources or destinations */
	*X = info;
	return (error ? -error : n_items);
}

int GMTMEX_post_process (void *API, struct GMTMEX *X, int n_items, mxArray *plhs[])
{	/* Get the data from GMT output items into the Matlab arrays */
	int item, k;
	unsigned int row, col;
	uint64_t gmt_ij;
	float *f = NULL;
	double *d = NULL;
	struct GMT_MATRIX *M = NULL;
	
	for (item = 0; item < n_items; item++) {
		if ((M = GMT_Retrieve_Data (API, X[item].ID)) == NULL) mexErrMsgTxt ("Error retrieving matrix from GMT\n");
		k = X[item].lhs_index;	/* Short-hand for index into plhs[] */
		switch (X[item].type) {
			case GMT_IS_GRID:	/* Return grids via float (mxSINGLE_CLASS) matrix */
				plhs[k] = mxCreateNumericMatrix (M->n_rows, M->n_columns, mxSINGLE_CLASS, mxREAL);
				f = mxGetData (plhs[k]);
				/* Load the real grd array into a double matlab array by transposing 
			           from unpadded GMT grd format to unpadded matlab format */
				for (gmt_ij = row = 0; row < M->n_rows; row++) for (col = 0; col < M->n_columns; col++, gmt_ij++)
					f[MEXM_IJ(M,row,col)] = M->data.f4[gmt_ij];
				break;
			case GMT_IS_DATASET:	/* Return tables with double (mxDOUBLE_CLASS) matrix */
				plhs[k] = mxCreateNumericMatrix (M->n_rows, M->n_columns, mxDOUBLE_CLASS, mxREAL);
				d = mxGetData (plhs[k]);
				/* Load the real data matrix into a double matlab array copying columns */
				for (col = 0; col < M->n_columns; col++)
					memcpy (&d[M->n_rows*col], M->data.f8, M->n_rows * sizeof(double));
				break;
		}
	}
	return(0);	/* Maybe we should turn this function to void but the gmt5 wraper expects a return value */
}
