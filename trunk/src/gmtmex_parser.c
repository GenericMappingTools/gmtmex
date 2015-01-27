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
#define GMTMEX_LIB
#include "gmtmex.h"
#include "gmtmex_modules.h"
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

/* We have one strsep in gmt.dll but it returns a temporary that is not visible here
   because of the cross DLLs shit. Se we use a local copy, named strsep_ */
#ifdef _MSC_VER
#	define STRSEP strsep_
#else
#	define STRSEP strsep
#endif

#ifndef rint
	#define rint(x) (floor((x)+0.5f)) //does not work reliable.
#endif

#if defined(WIN32) && !defined(lrint)
#	define lrint (int64_t)rint
#endif

#ifdef NO_MEX
/* For testing, no data are actually touched and Matlab/Octave is not linked */
int ID_lhs = 100000;
int ID_rhs = 0;
#define mxArray void
#define mexErrMsgTxt(txt) fprintf (stderr, txt)
#endif

/* New parser for all GMT mex modules based on design discussed by PW and JL on Mon, 2/21/11 */
/* Wherever we say "Matlab" we mean "Matlab of Octave" */

/* For the Mex interface we will wish to pass either filenames or matrices via GMT command options.
 * We select a Matlab matrix by suppying $ as the file name.  The parser will then find these $
 * arguments and replace them with references to matrices via the GMT API mechanisms.
 * This requires us to know which options in a module may accept a file name.  As an example,
 * consider surface whose -L option may take a grid.  To pass a Matlab/Octave grid already in memory
 * we would use -L$ and give the grid as an argument to the module, e.g.,
 *    Z = gmt ('surface -R0/50/0/50 -I1 -V xyzfile -L$', lowmatrix);
 * For each option that may take a file we need to know what kind of file and if this is input or output.
 * We encode this in a 3-character word, where
 *	1. The first char is the option flag (e.g., L for -L)
 *	2. The second is the data type (P|L|D|G|C|T)
 *	3. The third is I(n) or O(out)
 * E.g., the surface example would have the word LGI.  The data types P|L|D|G|C|T stand for
 * P(olygons), L(ines), D(point data), G(rid), C(PT file), T(ext table). [We originally only had
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

#define GMT_IS_PS		99	/* Use for PS output; use GMT_IS_GRID or GMT_IS_DATASET for data */


#ifdef GMT_MATLAB
/* Macros for getting the Matlab ij that correspond to (row,col) [no pad involved] */
/* This one operates on GMT_MATRIX */
#define MEXM_IJ(M,row,col) ((col)*M->n_rows + (row))
/* And this on GMT_GRID */
#define MEXG_IJ(M,row,col) ((col)*M->header->ny + M->header->ny - (row) - 1)
#else
/* Macros for getting the Octave ij that correspond to (row,col) [no pad involved] */
/* This one operates on GMT_MATRIX */
#define MEXM_IJ(M,row,col) ((row)*M->n_columns + (col))
/* And this on GMT_GRID */
#define MEXG_IJ(M,row,col) ((row)*M->header->nx + (col))
#endif

#ifdef _MSC_VER
char *strsep_ (char **stringp, const char *delim) {
	char *start = *stringp;
	char *ptr;

	if (start == NULL) return NULL;

	/* Optimize the case of no delimiters.  */
	if (delim[0] == '\0') {
		*stringp = NULL;
		return start;
	}

	/* Optimize the case of one delimiter.  */
	if (delim[1] == '\0')
		ptr = strchr (start, delim[0]);
	else
		/* The general case.  */
		ptr = strpbrk (start, delim);

	if (ptr == NULL) {
		*stringp = NULL;
		return start;
	}

	*ptr = '\0';
	*stringp = ptr + 1;

	return start;
}
#endif

int GMTMEX_find_module (void *API, char *module)
{	/* Just search for module and return entry in keys array.  Only modules listed in mexproginfo.txt are used */
	char gmt_module[GMT_STR16] = {"gmt"};
	int k, id = -1;
	for (k = 0; id == -1 && k < N_GMT_MODULES; k++)
		if (!strcmp (module, module_name[k]))
			id = k;
	if (id == -1) {	/* Not found in the known list, try prepending gmt to the module name (i.e. gmt + get = gmtget) */
		strcat (gmt_module, module);
		for (k = 0; id == -1 && k < N_GMT_MODULES; k++)
			if (!strcmp (gmt_module, module_name[k]))
				id = k;
		if (id == -1) return (-1);	/* Not found in the known list */
	}
	/* OK, found in the list - now call it and see if it is actually available */
	if ((k = GMT_Call_Module (API, module_name[id], GMT_MODULE_EXIST, NULL)) == GMT_NOERROR)	/* Found and accessible */
		return (id);
	return (-1);	/* Not found in any shared libraries */
}

#ifdef NO_MEX
#define mxstrdup(s) strdup(s)
#else
char *mxstrdup (const char *s) {
	/* A strdup replacement to be used in Mexs to avoid memory leaks since the Matlab
	   memory management will take care to free the memory allocated by this function */
	char *d = mxMalloc (strlen (s) + 1);
	if (d == NULL) return NULL;
	strcpy (d,s);
	return d;
}

int GMTMEX_print_func (FILE *fp, const char *message)
{
	/* Replacement for GMT's gmt_print_func.  It is being used indirectly via
	 * API->print_func.  Purpose of this is to allow Matlab (which cannot use
	 * printf) to reset API->print_func to this function via GMT_Create_Session. */

	mexPrintf (message);
	return 0;
}
#endif

int gmtmex_find_option (char option, char *key[], int n_keys) {
	/* gmtmex_find_option determines if the given option is among the special options listed
           in the key array that might take $ as filename */
	int pos = -1, k;
	for (k = 0; pos == -1 && k < n_keys; k++)
		if (key[k][0] == option) pos = k;
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

unsigned int gmtmex_get_key_pos (char *key[], unsigned int n_keys, struct GMT_OPTION *head, int def[2][2])
{	/* Must determine if default input and output have been set via program options or if they should be added explicitly.
 	 * As an example, consider the GMT command grdfilter in.nc -Fg200k -Gfilt.nc.  In Matlab this might be
	 * filt = gmt ('grdfilter $ -Fg200k -G$', in);
	 * However, it is more natural not to specify the lame -G$, i.e.
	 * filt = gmt ('grdfilter $ -Fg200k', in);
	 * or even the other lame $, e.g.
	 * filt = gmt ('grdfilter -Fg200k', in);
	 * In that case we need to know that -G is the default way to specify the output grid and if -G is not given we
	 * must associate -G with the first left-hand-side item (here filt).
	 */
	int pos, dir, flavor, PS = 0;
	struct GMT_OPTION *opt = NULL;
	def[GMT_IN][0] = GMT_MEX_IMPLICIT;	/* Initialize to setting the i/o implicitly for filenames */
	def[GMT_OUT][0] = GMT_MEX_NONE;		/* Initialize to setting the i/o implicitly for filenames */
	def[GMT_IN][1] = def[GMT_OUT][1] = GMT_MEX_NONE;	/* For options with mising filenames they are NONE unless set  */

	/* Loop over the module options to see if inputs and outputs are set explicitly or implicitly */
	for (opt = head; opt; opt = opt->next) {
		pos = gmtmex_find_option (opt->option, key, n_keys);	/* First see if this option is one that might take $ */
		if (pos == -1) continue;	/* No, it was some other harmless option, e.g., -J, -O ,etc. */
		//flavor = (opt->option == '<' || opt->option == '>') ? 0 : 1;	/* Filename or option with filename ? */
		flavor = (opt->option == '<') ? 0 : 1;	/* Filename or option with filename ? */
		dir = (key[pos][2] == 'I') ? GMT_IN : GMT_OUT;	/* Input of output ? */
		if (flavor == 0)	/* File name was given on command line */
			def[dir][flavor] = GMT_MEX_EXPLICIT;
		else			/* Command option; e.g., here we have -G<file>, -G$, or -G [the last two means implicit] */
			def[dir][flavor] = (opt->arg[0] == '\0' || opt->arg[0] == '$') ? GMT_MEX_IMPLICIT : GMT_MEX_EXPLICIT;	/* The option provided no file name (or gave $) so it is implicit */
	}
	/* Here, if def[] == GMT_MEX_IMPLICIT (the default in/out option was NOT given),
           then we want to return the corresponding entry in key */
	for (pos = 0; pos < n_keys; pos++) {	/* For all module options that might take a file */
		//flavor = (key[pos][0] == '<' || key[pos][0] == '>') ? 0 : 1;
		flavor = (key[pos][0] == '<') ? 0 : 1;
		if ((key[pos][2] == 'I' || key[pos][2] == 'i') && key[pos][0] == '-')
			/* This program takes no input (e.g., psbasemap, pscoast) */
			def[GMT_IN][0] = def[GMT_IN][1]  = GMT_MEX_NONE;
		else if (key[pos][2] == 'I' && def[GMT_IN][flavor] == GMT_MEX_IMPLICIT)
			/* Must add implicit input; use def to determine option,type */
			def[GMT_IN][flavor] = pos;
		else if ((key[pos][2] == 'O' || key[pos][2] == 'o') && key[pos][0] == '-')
			/* This program produces no output */
			def[GMT_OUT][0] = def[GMT_OUT][1] = GMT_MEX_NONE;
		else if (key[pos][2] == 'O' && def[GMT_OUT][flavor] == GMT_MEX_IMPLICIT)
			/* Must add implicit output; use def to determine option,type */
			def[GMT_OUT][flavor] = pos;
		else if (key[pos][2] == 'O' && def[GMT_OUT][flavor] == GMT_MEX_NONE && flavor == 1)
			/* Must add mising output option; use def to determine option,type */
			def[GMT_OUT][flavor] = pos;
		if ((key[pos][2] == 'O' || key[pos][2] == 'o') && key[pos][1] == 'X' && key[pos][0] == '-') PS = 1;	/* This program produces PostScript */
	}
	return (PS);
}

int gmtmex_get_arg_dir (char option, char *key[], int n_keys, int *data_type, int *geometry)
{	/* key[] is an array with options of the current program that read/write data */
	int item;

	/* 1. First determine if option is one of the choices in key */

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
	return ((key[item][2] == 'i') ? GMT_IN : GMT_OUT);	/* Return the direction of i/o */
}

char **make_char_array (char *string, unsigned int *n_items, char type)
{	/* Turn the comma-separated list of 3-char codes into an array of such codes.
 	 * In the process, replace any ?-types with the selected type. */
	size_t len, k, n;
	char **s = NULL;
	char *next, *tmp;

	if (!string) return NULL;	/* Got NULL, just give up */
	len = strlen (string);
	if (len == 0) return NULL;	/* Got no characters, give up */
	tmp = strdup (string);		/* Get a working copy of string */
	/* Replace unknown types in tmp with selected type */
	if (type)
		for (k = 0; k < strlen (tmp); k++)
			if (tmp[k] == '?') tmp[k] = type;
	/* Count the number of items */
	for (k = n = 0; k < len; k++)
		if (tmp[k] == ',') n++;
	n++;	/* Since one less comma than items */
	/* Allocate and populate the character array, then return it and n_items */
	s = (char **) calloc (n, sizeof (char *));
	k = 0;
	while ((next = STRSEP(&tmp, ",")) != NULL) {
		s[k++] = strdup (next);
	}
	*n_items = (unsigned int)n;
	free ((void *)tmp);
	return s;
}

#ifdef NO_MEX
int GMTMEX_Register_IO (void *API, unsigned int data_type, unsigned int geometry, unsigned int direction, const mxArray *ptr)
{
	if (direction == GMT_IN)
		return (ID_rhs++);	/* Fake IDs */
	else
		return (ID_lhs++);	/* Fake IDs */
}
#else
struct GMT_GRID *GMTMEX_grid_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Grid container to hold a GMT grid.
 	 * If direction is GMT_IN then we are given a Matlab grid and can determine its size, etc.
	 * If direction is GMT_OUT then we allocate an empty GMT grid that we will pass
	 * off as a destination by adding the GMT_VIA_OUTPUT to the mode. */
	struct GMT_GRID *G = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		mxArray *mx_ptr = NULL;
		double *inc = NULL, *range = NULL, *reg = NULL;
		unsigned int registration;

		if (mxIsEmpty (ptr))
			mexErrMsgTxt ("The input that was supposed to contain the Grid, is empty\n");
		if (!mxIsStruct (ptr))
			mexErrMsgTxt ("Expected a Grid for input\n");
		mx_ptr = mxGetField (ptr, 0, "inc");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("Could not find inc array with Grid increments\n");
		inc = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "range");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("Could not find range array for Grid range\n");
		range = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "registration");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("Could not find registration array for Grid registration\n");
		reg = mxGetData (mx_ptr);
		registration = (unsigned int)lrint (reg[0]);
		if ((G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY,
					NULL, range, inc, registration, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_grid_init: Failure to alloc GMT source matrix for input\n");
		mx_ptr = mxGetField (ptr, 0, "z");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("Could not find data array for Grid\n");

		G->data = mxGetData (mx_ptr);
		G->alloc_mode = GMT_ALLOCATED_EXTERNALLY;	/* Since array was allocated by Matlab */
		GMT_Report (API, GMT_MSG_DEBUG, " Allocate GMT Grid %lx in gmtmex_parser\n", (long)G);
		GMT_Report (API, GMT_MSG_DEBUG, " Registered GMT Grid array %lx via memory reference from Matlab\n", (long)G->data);
	}
	else {	/* Just allocate an empty container to hold the output grid, and pass GMT_VIA_OUTPUT */
		if ((G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY |
                         GMT_VIA_OUTPUT, NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_grid_init: Failure to alloc GMT blank grid container for holding output grid\n");
	}
	return (G);
}

struct GMT_MATRIX *GMTMEX_matrix_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Matrix container and associate it with a data matrix.
 	 * Note that in GMT these will be considered DATASETs via GMT_MATRIX.
 	 * If direction is GMT_IN then we are given a Matlab matrix and can determine size etc.
	 * If output then we dont know size but we can specify type */
	uint64_t dim[3] = {0, 0, 0}, *this_dim = NULL;
	unsigned int mode = 0;
	struct GMT_MATRIX *M = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known */
		if (!mxIsNumeric (ptr)) mexErrMsgTxt ("Expected a Matrix for input\n");
		dim[0] = mxGetN (ptr);
		dim[1] = mxGetM (ptr);
		this_dim = dim;
	}
	else	/* There are no dimensions yet, as we are getting data as output */
		mode = GMT_VIA_OUTPUT;	/* And this_dim == NULL */
	if ((M = GMT_Create_Data (API, GMT_IS_MATRIX, GMT_IS_PLP, mode, this_dim, NULL, NULL, 0, 0, NULL)) == NULL)
		mexErrMsgTxt ("GMTMEX_matrix_init: Failure to alloc GMT source matrix\n");

	GMT_Report (API, GMT_MSG_DEBUG, " Allocate GMT Matrix %lx in gmtmex_parser\n", (long)M);
	M->n_rows    = dim[1];
	M->n_columns = dim[0];
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
		/* Data from Matlab is in col format and from Octave in row format */
#ifdef GMT_MATLAB
		M->dim = M->n_rows;
#else
		M->dim = M->n_columns;
#endif
		M->alloc_mode = GMT_ALLOCATED_EXTERNALLY;	/* Since matrix was allocated by Matlab */
		M->shape = MEX_COL_ORDER;			/* Either col or row order, depending on Matlab/Octave */
	}
	else {	/* On output we produce double precision */
		M->type = GMT_DOUBLE;
		/* Data from GMT must be in row format since we may not know n_rows until later! */
		M->shape = GMT_IS_ROW_FORMAT;
	}

	return (M);
}

int GMTMEX_Register_IO (void *API, unsigned int data_type, unsigned int geometry, unsigned int direction, const mxArray *ptr)
{	/* Create the grid or matrix contains, register them, and return the ID */
	int ID = GMT_NOTSET;
	struct GMT_GRID *G = NULL;	/* Pointer to grid container */
	struct GMT_MATRIX *M = NULL;	/* Pointer to matrix container */

	switch (data_type) {
		case GMT_IS_GRID:
			/* Get an empty grid, and if input we and associate it with the Matlab grid pointer */
			G = GMTMEX_grid_init (API, direction, ptr);
			ID = GMT_Get_ID (API, GMT_IS_GRID, direction, G);
			GMT_Insert_Data (API, ID, G);
			break;
		case GMT_IS_DATASET:
			/* Get a matrix container, and if input and associate it with the Matlab pointer */
			M = GMTMEX_matrix_init (API, direction, ptr);
			ID = GMT_Get_ID (API, GMT_IS_DATASET, direction, M);
			break;
		default:
			mexErrMsgTxt ("GMTMEX_pre_process: Bad data type\n");
			break;
	}
	return (ID);
}
#endif

int GMTMEX_pre_process (void *API, const char *module, mxArray *plhs[], int nlhs, const mxArray *prhs[],
                        int nrhs, char *keys, struct GMT_OPTION **head, struct GMTMEX **X) {
	/* API controls all things within GMT.
	 * plhs (and nlhs) are the outputs specified on the left side of the equal sign in Matlab/Octave.
	 * prhs (and nrhs) are the inputs specified after the option string in the gmt call.
	 * keys is a comma-separated string with 3-char codes for current module i/o.
	 * opt is the linked list of GMT options passed in. X is a returned array of structures with
	 * information about registered resources going to/from GMT.
	 */

	int lr_pos[2] = {0, 0};	/* These position keeps track where we are in the L and R pointer arrays */
	int direction;		/* Either GMT_IN or GMT_OUT */
	int data_type;		/* Either GMT_IS_DATASET, GMT_IS_TEXTSET, GMT_IS_GRID, GMT_IS_CPT, GMT_IS_IMAGE */
	int geometry;		/* Either GMT_IS_NONE, GMT_IS_POINT, GMT_IS_LINE, GMT_IS_POLY, or GMT_IS_SURFACE */
	int given[2][2];	/* Either GMT_MEX_EXPLICIT or the item number in the keys array for a filename or an option with implicit arg */
	int ID, error;
	unsigned int k, flavor, n_keys = 0, PS, n_alloc = 8U, n_items = 0;
	char name[GMT_STR16];	/* Used to hold the GMT API embedded file name, e.g., @GMTAPI@-###### */
	char **key = NULL;
	char *text = NULL;
	char type = 0;
	void *ptr = NULL;
	struct GMT_OPTION *opt, *new_ptr;	/* Pointer to a GMT option structure */
	struct GMTMEX *info = NULL;
#ifndef NO_MEX
	struct GMT_GRID *G = NULL;		/* Pointer to grid container */
	struct GMT_MATRIX *M = NULL;		/* Pointer to matrix container */
#endif
	struct GMTAPI_CTRL *APIPI = NULL;
	APIPI = GMT_get_API_ptr(API);

	/* First, we check if this is either the read of write special module, which specifies what data type to deal with */
	if (!strcmp (module, "read") || !strcmp (module, "gmtread") || !strcmp (module, "write") || !strcmp (module, "gmtwrite")) {
		/* Special case: Must determine which data type we are dealing with */
		struct GMT_OPTION *t_ptr = NULL;
		if ((t_ptr = GMT_Find_Option (API, 'T', *head))) {	/* Found the -T<type> option */
			type = toupper (t_ptr->arg[0]);	/* Find type and replace ? in keys with this type in uppercase (DGCIT) in make_char_array below */
		}
		if (!strchr ("DGCIT", type)) {
			mexErrMsgTxt ("GMTMEX_pre_process: No or bad data type given to read|write\n");
		}
		if (!strstr ("write", module) && (t_ptr = GMT_Find_Option (API, GMT_OPT_INFILE, *head))) {
			/* Found a -<<file> option; this is actually the output file */
			t_ptr->option = GMT_OPT_OUTFILE;
		}
	}

	key = make_char_array (keys, &n_keys, type);		/* This is the array of keys for this module, e.g., "<DI,GGO,..." */
	info = malloc (n_alloc * sizeof (struct GMTMEX));	/* Structure to keep track of which output items we need to assign to Matlab */

	/* We wish to enable Matlab/Octave by "implicit" options.  These are options we will provide here in the code
	 * when the user did not specifically give them as part of the command.  For instance, if surface is called and the input data
	 * comes from a Matlab matrix, we may leave off any input name in the command, and then it is understood
	 * that we need to add a memory reference to the first Matlab matrix given as input. */

	PS = gmtmex_get_key_pos (key, n_keys, *head, given);	/* Determine if we must add the primary in and out arguments to the option list */
	/* Note: PS will be one if this module produces PostScript */
	for (direction = GMT_IN; direction <= GMT_OUT; direction++) {	/* Separately consider input and output */
		for (flavor = 0; flavor < 2; flavor++) {	/* 0 means filename input, 1 means option input */
			if (given[direction][flavor] == GMT_MEX_NONE) continue;		/* No source or destination required by this module */
			if (given[direction][flavor] == GMT_MEX_EXPLICIT) continue;	/* Source or destination was set explicitly in the command; skip */
			/* Here we must add the primary input or output from prhs[0] or plhs[0] */
			/* Get info about the data set */
			if (given[direction][flavor] < 0)
				mexErrMsgTxt("GMTMEX_pre_process: stoping here instead of crashing Matlab.\n\t\t'given[direction][flavor]' is negative\n");

			(void)gmtmex_get_arg_dir (key[given[direction][flavor]][0], key, n_keys, &data_type, &geometry);
			/* Pick the next left or right side Matlab array pointer */
			ptr = (direction == GMT_IN) ? (void *)prhs[lr_pos[GMT_IN]] : (void *)plhs[lr_pos[GMT_OUT]];
			/* Create and thus register this container */
			ID = GMTMEX_Register_IO (API, data_type, geometry, direction, ptr);
			if (ID == GMT_NOTSET)
				mexErrMsgTxt("GMTMEX_pre_process: Failure to register the resource\n");

			/* Keep a record or this container as a source or destination */
			if (n_items == n_alloc) info = realloc ((void *)info, (n_alloc += 8) * sizeof (struct GMTMEX));
			info[n_items].type = data_type;
			info[n_items].ID = ID;
			info[n_items].direction = direction;
			info[n_items].lhs_index = lr_pos[direction];
			n_items++;
			lr_pos[direction]++;		/* Advance position counter for next time */
			if (GMT_Encode_ID (API, name, ID) != GMT_NOERROR) 	/* Make filename with embedded object ID */
				mexErrMsgTxt ("GMTMEX_pre_process: Failure to encode string\n");

			if (flavor == 0) {	/* Must add a new option */
				/* Create the missing (implicit) GMT option */
				new_ptr = GMT_Make_Option (API, key[given[direction][0]][0], name);
				/* Append it to the option list */
				*head = GMT_Append_Option (API, new_ptr, *head);
			}
			else {	/* Must find the option and update it, or add it if not found */
				if ((new_ptr = GMT_Find_Option (API, key[given[direction][1]][0], *head)) == NULL) {
					/* Create the missing (implicit) GMT option */
					new_ptr = GMT_Make_Option (API, key[given[direction][1]][0], name);
					/* Append it to the option list */
					*head = GMT_Append_Option (API, new_ptr, *head);
				}
				else if (GMT_Update_Option (API, new_ptr, name)) {	/* Just update its arbument */
					mexErrMsgTxt ("GMTMEX_pre_process: Failure to update option argument\n");
				}
			}
		}
	}

	for (opt = *head; opt; opt = opt->next) 	/* Loop over the module options given */
		if (PS && opt->option == GMT_OPT_OUTFILE) PS++;	/* Count additional output options */

	/* Reallocate the information structure array or remove entirely if nothing given. */
	if (n_items && n_items < n_alloc) info = realloc ((void *)info, n_items * sizeof (struct GMTMEX));
	else if (n_items == 0) free ((void *)info);	/* No containers used */

	if (PS == 1)	/* No redirection of the PS to an actual file means an error */
		error = GMT_NOERROR;
	else if (PS > 2)	/* Too many output files for PS */
		error = 2;
	else
		error = GMT_NOERROR;
	/* Free up the temporary key array */
	for (k = 0; k < n_keys; k++) free ((void *)key[k]);
	free ((void *)key);

	/* Just checking that the options were properly */
	text = GMT_Create_Cmd (API, *head);
#ifdef NO_MEX
	GMT_Report (API, GMT_MSG_NORMAL, "Revised mex string for %s is \'%s\'\n", module, text);
#else
	GMT_Report (API, GMT_MSG_VERBOSE, "Args are now [%s]\n", text);
#endif
	GMT_Destroy_Cmd (API, &text);

	/* Here, a command line '-F200k -G $ -L -P' has been changed to '-F200k -G@GMTAPI@-000001 @GMTAPI@-000002 -L@GMTAPI@-000003 -P'
	 * where the @GMTAPI@-00000x are encodings to registered resources or destinations */

	/* Pass back the info array and the number of items */
	*X = info;
	return (error ? -error : n_items);
}

#ifdef NO_MEX
int GMTMEX_post_process (void *API, struct GMTMEX *X, int n_items, mxArray *plhs[]) {
	GMT_Report (API, GMT_MSG_VERBOSE, "Exit GMTMEX_post_process\n");
	return (0);
}
#else
int GMTMEX_post_process (void *API, struct GMTMEX *X, int n_items, mxArray *plhs[]) {
	/* Get the data from GMT output items into the corresponding Matlab struct or matrix */
	int item, k, n;
	unsigned int row, col;
	uint64_t gmt_ij, mex_ij;
	float  *f = NULL;
	double *d = NULL, *dptr = NULL, *G_x = NULL, *G_y = NULL, *x = NULL, *y = NULL;
	struct GMT_GRID *G = NULL;
	struct GMT_MATRIX *M = NULL;
	mxArray *mxGrd = NULL, *mx_x = NULL, *mx_y= NULL;
	mxArray *mxProjectionRef = NULL;
	mxArray *mxHeader = NULL, *mxtmp = NULL;
	mxArray *grid_struct = NULL;
	char    *fieldnames[22];	/* this array contains the names of the fields of the output grid structure. */

	memset ((void *)fieldnames, 0, 22*sizeof (char *));
	GMT_Report (API, GMT_MSG_VERBOSE, "Enter GMTMEX_post_process\n");

	for (item = 0; item < n_items; item++) {	/* Number of GMT container involved in the call */
		k = X[item].lhs_index;	/* Short-hand for index into plhs[] */
		switch (X[item].type) {
			case GMT_IS_GRID:	/* We read or wrote a GMT grid, examine further */
				if ((G = GMT_Retrieve_Data (API, X[item].ID)) == NULL)
					mexErrMsgTxt ("Error retrieving grid from GMT\n");
				if (X[item].direction == GMT_OUT) {	/* Here, GMT_OUT means "Return this info to Matlab" */
					/* Return grids via a float (mxSINGLE_CLASS) matrix in a struct */
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
					fieldnames[16] = mxstrdup ("x");
					fieldnames[17] = mxstrdup ("y");
					fieldnames[18] = mxstrdup ("z");
					fieldnames[19] = mxstrdup ("x_units");
					fieldnames[20] = mxstrdup ("y_units");
					fieldnames[21] = mxstrdup ("z_units");
					grid_struct = mxCreateStructMatrix (1, 1, 22, (const char **)fieldnames );

					mxtmp = mxCreateString (G->header->ProjRefPROJ4);
					mxSetField (grid_struct, 0, (const char *) "ProjectionRefPROJ4", mxtmp);

					mxtmp = mxCreateString (G->header->ProjRefWKT);
					mxSetField (grid_struct, 0, (const char *) "ProjectionRefWKT", mxtmp);

					mxHeader    = mxCreateNumericMatrix (1, 9, mxDOUBLE_CLASS, mxREAL);
					dptr = mxGetPr (mxHeader);
					for (n = 0; n < 4; n++) dptr[n] = G->header->wesn[n];
					dptr[4] = G->header->z_min;	dptr[5] = G->header->z_max;
					dptr[6] = G->header->registration;
					for (n = 0; n < 2; n++) dptr[n+7] = G->header->inc[n];
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

					mxtmp = mxCreateString ("float32");
					mxSetField (grid_struct, 0, (const char *) "DataType", mxtmp);

					mxtmp = mxCreateDoubleScalar ((double)G->header->n_bands);
					mxSetField (grid_struct, 0, (const char *) "LayerCount", mxtmp);

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
					for (gmt_ij = row = 0; row < G->header->ny; row++)
						for (col = 0; col < G->header->nx; col++, gmt_ij++)
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
					for (n = 0; n < G->header->ny; n++) y[G->header->ny-1-n] = G_y[n];
					if (GMT_Destroy_Data (API, &G_x))
						mexPrintf("Warning: Failure to delete G_x (x coordinate vector)\n");
					if (GMT_Destroy_Data (API, &G_y))
						mexPrintf("Warning: Failure to delete G_y (y coordinate vector)\n");
					mxSetField (grid_struct, 0, "x", mx_x);
					mxSetField (grid_struct, 0, "y", mx_y);
					plhs[k] = grid_struct;
				}
				/* Else it is a matlab grid that was passed into GMT as input */
				/* We always destroy G at this point, whether input or output.  The alloc_mode
				 * will prevent accidential freeing of any Matlab-allocated arrays */
				if (GMT_Destroy_Data (API, &G) != GMT_NOERROR)
					mexErrMsgTxt ("GMTMEX_post_process: Failed to destroy grid G used in the interface bewteen GMT and Matlab\n");
				break;

			case GMT_IS_DATASET:	/* Return tables with double (mxDOUBLE_CLASS) matrix */
				if (X[item].direction == GMT_OUT) {		/* Here, GMT_OUT means "Return this info to Matlab" */
					if ((M = GMT_Retrieve_Data(API, X[item].ID)) == NULL)
						mexErrMsgTxt ("Error retrieving matrix from GMT\n");

					/* Create a Matlab matrix to hold this GMT matrix/data table */
					plhs[k] = mxCreateNumericMatrix (M->n_rows, M->n_columns, mxDOUBLE_CLASS, mxREAL);
					d = mxGetData (plhs[k]);
					/* Duplicate the double data matrix into the matlab double array */
					if (M->shape == MEX_COL_ORDER)	/* Easy, just copy */
						memcpy (d, M->data.f8, M->n_rows * M->n_columns * sizeof (double));
					else {	/* Must transpose */
						for (gmt_ij = row = 0; row < M->n_rows; row++) {
							for (col = 0; col < M->n_columns; col++, gmt_ij++) {
								mex_ij = MEXM_IJ (M, row, col);
								d[mex_ij] = M->data.f8[gmt_ij];
							}
						}
					}
			//		/* Load the real data matrix into a double matlab array copying columns */
			//		for (col = 0; col < M->n_columns; col++)
			//			memcpy (&d[M->n_rows*col], M->data.f8, M->n_rows * sizeof(double));
				}
				/* Else we were passing Matlab data into GMT as data input */
				/* We always destroy M at this point, whether input or output.  The alloc_mode
				 * will prevent accidential freeing of any Matlab-allocated arrays */
				if (GMT_Destroy_Data (API, &M) != GMT_NOERROR)
					mexErrMsgTxt ("GMTMEX_post_process: Failed to destroy matrix M used in the interface bewteen GMT and Matlab\n");

				break;
		}
	}
	GMT_Report (API, GMT_MSG_VERBOSE, "Exit GMTMEX_post_process\n");
	return (0);	/* Maybe we should turn this function to void but the gmt5 wrapper expects a return value */
}
#endif
