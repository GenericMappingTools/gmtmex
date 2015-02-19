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

#if defined(WIN32)
#if !defined(lrint)
#	define lrint (int64_t)rint
#endif
#include <io.h>
#else
#include <unistd.h>
#endif
#ifndef F_OK
#	define F_OK 0
#endif
#ifdef NO_MEX
/* For testing, no data are actually touched and Matlab/Octave is not linked */
int ID_lhs = 100000;	/* We start output IDs at 100000 for fun */
int ID_rhs = 0;
#define mxArray void
#define mexErrMsgTxt(txt) fprintf (stderr, txt)
#endif

/* Indices into the keys triple code */
#define K_OPT	0
#define K_TYPE	1
#define K_DIR	2

enum MEX_dim {
	DIM_COL	= 0,	/* Holds the number of columns for vectors and x-nodes for matrix */
	DIM_ROW = 1};	/* Holds the number of rows for vectors and y-nodes for matrix */

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
 * We encode this in a 3-character word XYZ, explained below.  Note that each module may
 * need several comma-separated XYZ words and these are returned as one string via GMT_Get_Moduleinfo.
 *
 * X stands for the specific program option (e.g., L for -L, F for -F), or <,>
 *    for standard input,output (if reading tables) or command-line files (if reading grids).
 *    A hyphen (-) means there is no option for this item.
 * Y stands for data type (C = CPT, D = Dataset/Point, L = Dataset/Line,
 *    P = Dataset/Polygon, G = Grid, I = Image, T = Textset, X = PostScript, ? = type given via module option),
 * Z stands for primary inputs (I), primary output (O), secondary input (i) secondary output (o).
 *   Primary inputs and outputs need to be assigned, and if not explicitly given we will
 *   use the given left- and right-hand side arguments to supply input or accept output.
 *   Secondary inputs means they are only assigned if the option is given.
 *
 * E.g., the surface example would have the word LGI.  The data types P|L|D|G|C|T stand for
 * P(olygons), L(ines), D(point data), G(rid), C(PT file), T(ext table). [We originally only had
 * D for datasets but perhaps the geometry needs to be passed too (if known); hence the P|L|D char]
 * In addition, the only common option that might take a file is -R which may take a grid as input.
 * We check for that in addition to the module-specific info passed via the key variable.
 *
 * The actual reading/writing will occur in gmt_api.c via the standard GMT containers.
 * The packing up GMT grids into Matlab grid structs and vice versa happens after getting the
 * results out of the GMT API and before passing into back to Matlab.
 */

/* We will consolidate this code once everything is working.  Parts of this code (the things that
 * only depend on GMT functions) will be included in the API and only documented in the API developer
 * section while things that is tied to the external languate (Matlab, Python, etc) will remain here.
 * We flag sections either as GMT_ONLY or EXTERNAL below for now. */

#define GMT_FILE_NONE		0
#define GMT_FILE_EXPLICIT	1
#define GMT_FILE_IMPLICIT	2

#define GMT_IS_PS		99	/* Use for PS output; use GMT_IS_GRID or GMT_IS_DATASET for data */

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

#define N_MEX_FIELDNAMES	22

void * GMTMEX_Get_Grid (void *API, struct GMT_GRID *G)
{	/* Hook this grid into the k'th output item */

	int item, n;
	unsigned int row, col;
	uint64_t gmt_ij, mex_ij;
	float  *f = NULL;
	double *d = NULL, *dptr = NULL, *G_x = NULL, *G_y = NULL, *x = NULL, *y = NULL;
	mxArray *mxGrd = NULL, *mx_x = NULL, *mx_y= NULL;
	mxArray *mxProjectionRef = NULL;
	mxArray *mxHeader = NULL, *mxtmp = NULL;
	mxArray *grid_struct = NULL;
	char    *fieldnames[N_MEX_FIELDNAMES];	/* this array contains the names of the fields of the output grid structure. */

	if (!G->data)
		mexErrMsgTxt ("GMTMEX_Get_Grid: programming error, output matrix G is empty\n");

	memset (fieldnames, 0, N_MEX_FIELDNAMES*sizeof (char *));
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
	grid_struct = mxCreateStructMatrix (1, 1, N_MEX_FIELDNAMES, (const char **)fieldnames );

	mxtmp = mxCreateString (G->header->ProjRefPROJ4);
	mxSetField (grid_struct, 0, (const char *) "ProjectionRefPROJ4", mxtmp);

	mxtmp = mxCreateString (G->header->ProjRefWKT);
	mxSetField (grid_struct, 0, (const char *) "ProjectionRefWKT", mxtmp);

	mxHeader = mxCreateNumericMatrix (1, 9, mxDOUBLE_CLASS, mxREAL);
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

	mxGrd = mxCreateNumericMatrix (G->header->ny, G->header->nx, mxSINGLE_CLASS, mxREAL);
	f = mxGetData (mxGrd);
	/* Load the real grd array into a double matlab array by transposing
           from unpadded GMT grd format to unpadded matlab format */
	for (gmt_ij = row = 0; row < G->header->ny; row++)
		for (col = 0; col < G->header->nx; col++, gmt_ij++)
			f[MEXG_IJ(G,row,col)] = G->data[gmt_ij];
	mxSetField (grid_struct, 0, "z", mxGrd);

	/* Also return the convenient x and y arrays */
	G_x = GMT_Get_Coord (API, GMT_IS_GRID, GMT_X, G);	/* Get array of x coordinates */
	G_y = GMT_Get_Coord (API, GMT_IS_GRID, GMT_Y, G);	/* Get array of y coordinates */
	mx_x = mxCreateNumericMatrix (1, G->header->nx, mxDOUBLE_CLASS, mxREAL);
	mx_y = mxCreateNumericMatrix (1, G->header->ny, mxDOUBLE_CLASS, mxREAL);
	x = mxGetData (mx_x);
	y = mxGetData (mx_y);
	memcpy (x, G_x, G->header->nx * sizeof (double));
	for (n = 0; n < G->header->ny; n++) y[G->header->ny-1-n] = G_y[n];	/* Must reverse the y-array */
	if (GMT_Destroy_Data (API, &G_x))
		mexPrintf("Warning: Failure to delete G_x (x coordinate vector)\n");
	if (GMT_Destroy_Data (API, &G_y))
		mexPrintf("Warning: Failure to delete G_y (y coordinate vector)\n");
	mxSetField (grid_struct, 0, "x", mx_x);
	mxSetField (grid_struct, 0, "y", mx_y);
	return (grid_struct);
}

void * GMTMEX_Get_Table (void *API, struct GMT_MATRIX *M) {	/* Hook this table into the k'th output item */
	unsigned int row, col;
	uint64_t gmt_ij, mex_ij;
	mxArray *P = mxCreateNumericMatrix (M->n_rows, M->n_columns, mxDOUBLE_CLASS, mxREAL);
	double *d = mxGetData (P);
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
	return (P);
}

void * GMTMEX_Get_Text (void *API, struct GMT_MATRIX *M) {
	GMT_Report (API, GMT_MSG_NORMAL, "INTERNAL ERROR: Handling of TEXTSET not implemented yet\n");
	return NULL;
}

void * GMTMEX_Get_CPT (void *API, struct GMT_MATRIX *M) {
	GMT_Report (API, GMT_MSG_NORMAL, "INTERNAL ERROR: Handling of CPT not implemented yet\n");
	return NULL;
}

void * GMTMEX_Get_Image (void *API, struct GMT_MATRIX *M) {
	GMT_Report (API, GMT_MSG_NORMAL, "INTERNAL ERROR: Handling of IMAGE not implemented yet\n");
	return NULL;
}
#endif

#ifdef NO_MEX
void * GMTMEX_Register_IO (void *API, unsigned int family, unsigned int geometry, unsigned int direction, const mxArray *ptr, int *ID)
{	/* For testing we do not hook anything up, just increment the fake IDs */
	if (direction == GMT_IN) {
		ID_rhs++;	/* Fake IDs */
		*ID = ID_rhs;
	}
	else {
		ID_lhs++;	/* Fake IDs */
		*ID = ID_lhs;
	}
	return (NULL);
}
#else
struct GMT_GRID *GMTMEX_grid_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Grid container to hold a GMT grid.
 	 * If direction is GMT_IN then we are given a Matlab grid and can determine its size, etc.
	 * If direction is GMT_OUT then we allocate an empty GMT grid as a destination. */
	unsigned int row, col;
	uint64_t gmt_ij;
	struct GMT_GRID *G = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		mxArray *mx_ptr = NULL;
		double *inc = NULL, *range = NULL, *reg = NULL, *MinMax = NULL;
		float *f = NULL;
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
		if ((G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL,
		                          NULL, range, inc, registration, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_grid_init: Failure to alloc GMT source matrix for input\n");

		mx_ptr = mxGetField (ptr, 0, "MinMax");
		if (mx_ptr != NULL) {	/* Because we sent a NULL instead of the data array, z_min, z_max are not known. Use those from ptr */
			MinMax = mxGetData (mx_ptr);
			G->header->z_min = MinMax[0];
			G->header->z_max = MinMax[1];
		}

		mx_ptr = mxGetField (ptr, 0, "z");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("Could not find data array for Grid\n");

		f = mxGetData (mx_ptr);
		for (gmt_ij = row = 0; row < G->header->ny; row++)
			for (col = 0; col < G->header->nx; col++, gmt_ij++)
				G->data [gmt_ij] = f[MEXG_IJ(G,row,col)];
		GMT_Report (API, GMT_MSG_DEBUG, " Allocate GMT Grid %lx in gmtmex_parser\n", (long)G);
		GMT_Report (API, GMT_MSG_DEBUG, " Registered GMT Grid array %lx via memory reference from Matlab\n", (long)G->data);
	}
	else {	/* Just allocate an empty container to hold an output grid (signal this by passing NULLs) */
		if ((G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY,
                        NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_grid_init: Failure to alloc GMT blank grid container for holding output grid\n");
	}
	return (G);
}

struct GMT_MATRIX *GMTMEX_matrix_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Matrix container and associate it with a data matrix.
 	 * Note that in GMT these will be considered DATASETs via GMT_MATRIX.
 	 * If direction is GMT_IN then we are given a Matlab matrix and can determine size, etc.
	 * If output then we dont know size but we can specify type */
	uint64_t dim[3] = {0, 0, 0}, *this_dim = NULL;
	struct GMT_MATRIX *M = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known, extract them and set this_dim pointer */
		if (!mxIsNumeric (ptr)) mexErrMsgTxt ("Expected a Matrix for input\n");
		dim[DIM_ROW] = mxGetM (ptr);	/* Number of rows */
		dim[DIM_COL] = mxGetN (ptr);	/* Number of columns */
		this_dim = dim;
	}
	/* Else there are no dimensions yet, this_dim is NULL and we are getting an empty container for output */
	if ((M = GMT_Create_Data (API, GMT_IS_MATRIX, GMT_IS_PLP, 0, this_dim, NULL, NULL, 0, 0, NULL)) == NULL)
		mexErrMsgTxt ("GMTMEX_matrix_init: Failure to alloc GMT source matrix\n");

	GMT_Report (API, GMT_MSG_DEBUG, " Allocate GMT Matrix %lx in gmtmex_parser\n", (long)M);
	M->n_rows    = dim[DIM_ROW];
	M->n_columns = dim[DIM_COL];
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
		/* Data from Matlab and Octave(mex) is in col format and data from Octave(oct) is in row format */
#ifdef GMT_OCTOCT
		M->dim = M->n_columns;
#else
		M->dim = M->n_rows;
#endif

		M->alloc_mode = GMT_ALLOC_EXTERNALLY;	/* Since matrix was allocated by Matlab/Octave */
		M->shape = MEX_COL_ORDER;		/* Either col or row order, depending on Matlab/Octave setting in gmtmex.h */
	}
	else {	/* On output we produce double precision */
		M->type = GMT_DOUBLE;
		/* Data from GMT must be in row format since we may not know n_rows until later! */
		M->shape = GMT_IS_ROW_FORMAT;
	}

	return (M);
}

void * GMTMEX_Register_IO (void *API, unsigned int family, unsigned int geometry, unsigned int direction, const mxArray *ptr, int *ID)
{	/* Create the grid or matrix container, register it, and return the ID */
	struct GMT_GRID *G = NULL;	/* Pointer to grid container */
	struct GMT_MATRIX *M = NULL;	/* Pointer to matrix container */
	void *obj = NULL;		/* Pointer to the container we created */
	*ID = GMT_NOTSET;

	switch (family) {
		case GMT_IS_GRID:
			/* Get an empty grid, and if input we and associate it with the Matlab grid pointer */
			G = GMTMEX_grid_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_GRID, direction, G);
			GMT_Insert_Data (API, *ID, G);
			obj = G;
			GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Got Grid G with ID %d\n", *ID);
			break;
		case GMT_IS_DATASET:
			/* Get a matrix container, and if input and associate it with the Matlab pointer */
			M = GMTMEX_matrix_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_DATASET, direction, M);
			GMT_Insert_Data (API, *ID, M);
			obj = M;
			GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Got Matrix M with ID %d\n", *ID);
			break;
		default:
			GMT_Report (API, GMT_MSG_LONG_VERBOSE, "GMTMEX_Register_IO: Bad data type (%d)\n", family);
			break;
	}
	return (obj);
}
#endif

/* These functions require GMT only and may go into gmt_api.c eventually.
 * Some may become part of the API:
 *   int GMT_Get_Info 		Parse the options and the module keys and build info array
 *   void GMT_Expand_Option	Replace a marker ($) with a memory filename
 * The others are called from GMT_Get_Info only.
 */ 

int GMTAPI_key_to_family (char *key, int *family, int *geometry)
{
	/* Assign direction, family, and geometry based on key */

	switch (key[K_TYPE]) {	/* 2nd char contains the data type code */
		case 'G':
			*family = GMT_IS_GRID;
			*geometry = GMT_IS_SURFACE;
			break;
		case 'P':
			*family = GMT_IS_DATASET;
			*geometry = GMT_IS_POLY;
			break;
		case 'L':
			*family = GMT_IS_DATASET;
			*geometry = GMT_IS_LINE;
			break;
		case 'D':
			*family = GMT_IS_DATASET;
			*geometry = GMT_IS_POINT;
			break;
		case 'C':
			*family = GMT_IS_CPT;
			*geometry = GMT_IS_NONE;
			break;
		case 'T':
			*family = GMT_IS_TEXTSET;
			*geometry = GMT_IS_NONE;
			break;
		case 'I':
			*family = GMT_IS_IMAGE;
			*geometry = GMT_IS_SURFACE;
			break;
		case 'X':
			*family = GMT_IS_PS;
			*geometry = GMT_IS_NONE;
			break;
		default:
			return -2;
			break;
	}

	/* Third key character contains the in/out code */
	return ((key[K_DIR] == 'i' || key[K_DIR] == 'I') ? GMT_IN : GMT_OUT);	/* Return the direction of the i/o */
}

char **GMTAPI_process_keys (void *API, const char *string, char type, unsigned int *n_items, unsigned int *PS)
{	/* Turn the comma-separated list of 3-char codes into an array of such codes.
 	 * In the process, replace any ?-types with the selected type if type is not 0. */
	size_t len, k, n;
	char **s = NULL, *next = NULL, *tmp = NULL;
	*PS = 0;	/* No PostScript output indicated so far */

	if (!string) return NULL;	/* Got NULL, just give up */
	len = strlen (string);		/* Get the length of this item */
	if (len == 0) return NULL;	/* Got no characters, give up */
	tmp = strdup (string);		/* Get a working copy of string */
	/* Replace unknown types (?) in tmp with selected type give by variable "type" */
	if (type)	/* Got a nonzero type */
		for (k = 0; k < strlen (tmp); k++)
			if (tmp[k] == '?') tmp[k] = type;
	/* Count the number of items */
	for (k = 0, n = 1; k < len; k++)
		if (tmp[k] == ',') n++;
	/* Allocate and populate the character array, then return it and n_items */
	s = (char **) calloc (n, sizeof (char *));
	k = 0;
	while ((next = STRSEP(&tmp, ",")) != NULL) {
		s[k++] = strdup (next);
		if (strlen (next) != 3)
			GMT_Report (API, GMT_MSG_NORMAL, "INTERNAL ERROR: key %s does not contain exactly 3 characters\n", next);
		if (!strcmp (next, "-Xo")) (*PS)++;	/* Found key for PostScript output */
	}
	*n_items = (unsigned int)n;
	free ((void *)tmp);
	return s;
}

int GMTAPI_get_key (char option, char *keys[], int n_keys)
{	/* Return the position in the keys array that matches this option, or -1 if not found */
	int k;
	for (k = 0; k < n_keys; k++) if (keys[k][K_OPT] == option) return (k);
	return (-1);
}

int GMT_Get_Info (void *API, char *module, char marker, struct GMT_OPTION **head, struct GMT_INFO **X) {
	/* This function determines which input sources and output destinations are required.
	 * These are the arguments:
	 *   API controls all things within GMT.
	 *   module is the name of the GMT module. While an input arg it may grow if a prefix of "gmt" is prepended.
	 *   marker is the character that represents a resource, typically $
	 *   head is the linked list of GMT options passed for this module.  We may hook on 1-2 additional options.
	 *   X is a returned array of structures with information about registered resources going to/from GMT.
	 *   Number of structures is returned by the function.
	 */

	unsigned int n_keys, direction, PS, kind;
	unsigned int n_items = 0, n_explicit = 0, n_implicit = 0, output_pos = 0, explicit_pos = 0, implicit_pos = 0;
	int family;		/* -1, or one of GMT_IS_DATASET, GMT_IS_TEXTSET, GMT_IS_GRID, GMT_IS_CPT, GMT_IS_IMAGE */
	int geometry;		/* -1, or one of GMT_IS_NONE, GMT_IS_POINT, GMT_IS_LINE, GMT_IS_POLY, GMT_IS_SURFACE */
	int k;
	size_t n_alloc, len;
	const char *keys = NULL;	/* This module's option keys */
	char **key = NULL;		/* Array of items in keys */
	char *text = NULL, *LR[2] = {"rhs", "lhs"}, *S[2] = {" IN", "OUT"}, txt[16] = {""};
	char type = 0;
	struct GMT_OPTION *opt = NULL, *new_ptr = NULL;	/* Pointer to a GMT option structure */
	struct GMT_INFO *info = NULL;	/* Our return array of n_items info structures */

	/* 0. Get the keys for the module, possibly prepend "gmt" to module if required, or list modules and return if unknown module */
	if ((keys = GMT_Get_Moduleinfo (API, module)) == NULL) {	/* Gave an unknown module */
		GMT_Call_Module (API, NULL, GMT_MODULE_PURPOSE, NULL);	/* List available modules */
		return (-1);	/* Unknown module */
	}

	/* 1. Check if this is either the read of write special module, which specifies what data type to deal with via -T<type> */
	if (!strncmp (module, "gmtread", 7U) || !strncmp (module, "gmtwrite", 8U)) {
		/* Special case: Must determine which data type we are dealing with via -T<type> */
		if ((opt = GMT_Find_Option (API, 'T', *head)))	/* Found the -T<type> option */
			type = toupper (opt->arg[0]);	/* Find type and replace ? in keys with this type in uppercase (DGCIT) in GMTAPI_process_keys below */
		if (!strchr ("DGCIT", type)) {
			GMT_Report (API, GMT_MSG_NORMAL, "No or bad data type given to read|write (%c)\n", type);
			return (-2);
		}
		if (!strncmp (module, "gmtwrite", 8U) && (opt = GMT_Find_Option (API, GMT_OPT_INFILE, *head))) {
			/* Found a -<"file" option; this is actually the output file so we reset the option */
			opt->option = GMT_OPT_OUTFILE;
		}
	}

	/* 2. Get the option key array for this module, and determine if it produces PostScript output (PS == 1) */
	key = GMTAPI_process_keys (API, keys, type, &n_keys, &PS);	/* This is the array of keys for this module, e.g., "<DI,GGO,..." */

	/* 3. Count the module options and any input files referenced via marker, then allocate info struct array */
	for (opt = *head; opt; opt = opt->next) {
		if (strchr (opt->arg, marker)) n_explicit++;	/* Found an explicit dollar sign referring to an input matrix */
		if (PS && opt->option == GMT_OPT_OUTFILE) PS++;	/* Count given output options when PS will be produced. */
	}
	if (PS == 1) {	/* No redirection of the PS to an actual file means an error */
		GMT_Report (API, GMT_MSG_NORMAL, "No PostScript output file given\n");
		return (-3);
	}
	else if (PS > 2) {
		GMT_Report (API, GMT_MSG_NORMAL, "Can only specify one PostScript output file\n");
		return (-4);
	}
	n_alloc = n_explicit + n_keys;	/* Max number of registrations needed (may be just n_explicit) */
	info = calloc (n_alloc, sizeof (struct GMT_INFO));

	/* 4. Determine position of file args given as $ or via missing arg (proxy for input matrix) */
	/* Note: All implicit options must be given after all implicit matrices have been listed */
	for (opt = *head, implicit_pos = n_explicit; opt; opt = opt->next) {	/* Process options */
		k = GMTAPI_get_key (opt->option, key, n_keys);	/* If k >= 0 then this option is among those listed in the keys array */
		family = geometry = -1;	/* Not set yet */
		if (k >= 0) direction = GMTAPI_key_to_family (key[k], &family, &geometry);	/* Get dir, datatype, and geometry */
		
		if (strchr (opt->arg, marker)) {	/* Found an explicit dollar sign [these are always inputs] */
			direction = GMT_IN;
			/* Note sure about the OPT_INFILE test - should apply to all, no? */
			if (k >= 0 && key[k][K_OPT] == GMT_OPT_INFILE) key[k][K_DIR] = tolower (key[k][K_DIR]);	/* Make sure required I becomes i so we dont add it later */
			info[n_items].option = opt;
			info[n_items].family = family;
			info[n_items].geometry = geometry;
			info[n_items].direction = direction;
			info[n_items].pos = explicit_pos++;	/* Explicitly given arguments are the first given on the r.h.s. */
			kind = GMT_FILE_EXPLICIT;
			n_items++;
		}
		else if (k >= 0 && key[k][K_OPT] != GMT_OPT_INFILE && (len = strlen (opt->arg)) < 2) {	/* Got some option like -G or -Lu with further args */
			/* This is an implicit reference and we must add the missing item */
			info[n_items].option = opt;
			info[n_items].family = family;
			info[n_items].geometry = geometry;
			info[n_items].direction = direction;
			key[k][K_DIR] = tolower (key[k][K_DIR]);	/* Change to lowercase i or o since option was provided, albeit implicitly */
			info[n_items].pos = (direction == GMT_IN) ? implicit_pos++ : output_pos++;
			/* Excplicitly add the missing marker ($) to the option argument */
			sprintf (txt, "%s%c", opt->arg, marker);
			free (opt->arg);
			opt->arg = strdup (txt);
			n_implicit++;
			kind = GMT_FILE_EXPLICIT;
			n_items++;
		}
		else
			kind = GMT_FILE_NONE;	/* No file argument involved */
		if (kind == GMT_FILE_EXPLICIT)
			GMT_Report (API, GMT_MSG_LONG_VERBOSE, "%s: Option -%c%s includes an explicit reference to %s argument # %d\n", S[direction], opt->option, opt->arg, LR[direction], info[n_items].pos);
		else
			GMT_Report (API, GMT_MSG_LONG_VERBOSE, "---: Option -%c%s includes no special arguments\n", opt->option, opt->arg);
	}
	
	/* Done processing things that were explicitly given in the options.  Now look if module has required
	 * input or output items that we must add if not specified already */
	
	for (k = 0; k < n_keys; k++) {	/* Each set of keys specifies if the item is required via the 3rd key letter */
		if (isupper (key[k][K_DIR])) {	/* Required, and was not reset above */
			char str[2] = {0,0};
			str[0] = marker;
			direction = GMTAPI_key_to_family (key[k], &family, &geometry);
			new_ptr = GMT_Make_Option (API, key[k][K_OPT], str);	/* Create new option with filename "$" */
			/* Append the new option to the list */
			*head = GMT_Append_Option (API, new_ptr, *head);
			info[n_items].option = new_ptr;
			info[n_items].family = family;
			info[n_items].geometry = geometry;
			info[n_items].direction = direction;
			info[n_items].pos = (direction == GMT_IN) ? implicit_pos++ : output_pos++;
			GMT_Report (API, GMT_MSG_LONG_VERBOSE, "%s: Must add -%c%c as implicit reference to %s argument # %d\n",
				S[direction], key[k][K_OPT], marker, LR[direction], info[n_items].pos);
			n_items++;
		}
		free ((void *)key[k]);	/* Free up this key */
	}
	/* Free up the temporary key array */
	free ((void *)key);

	/* Reallocate the information structure array or remove entirely if nothing given. */
	if (n_items && n_items < n_alloc) info = realloc ((void *)info, n_items * sizeof (struct GMT_INFO));
	else if (n_items == 0) free ((void *)info);	/* No containers used */


	GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Found %d inputs and %d outputs\n", implicit_pos, output_pos);
	/* Just checking that the options were properly processed */
	text = GMT_Create_Cmd (API, *head);
#ifdef NO_MEX
	sprintf (revised_cmd, "\'%s %s\'", module, text);
#else
	GMT_Report (API, GMT_MSG_LONG_VERBOSE, "Revised command: %s\n", text);
#endif
	GMT_Destroy_Cmd (API, &text);	/* Only needed it for the above verbose */

	/* Pass back the info array and the number of items */
	*X = info;
	return (n_items);
}

void GMT_Expand_Option (void *API, struct GMT_OPTION *option, char marker, char *txt)
{	/* Replace marker with txt in the option argument */
	char buffer[BUFSIZ] = {""};
	size_t i = 0, o = 0;
	while (option->arg[i]) {
		if (option->arg[i] == marker) {
			strcat (&buffer[o], txt);
			o += strlen (txt);
		}
		else
			buffer[o++] = option->arg[i];
		i++;
	}
	free (option->arg);
	option->arg = strdup (buffer);
}
