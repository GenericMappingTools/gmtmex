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

#ifndef rint
	#define rint(x) (floor((x)+0.5f)) //does not work reliable.
#endif

#if defined(WIN32)
#define strdup _strdup
#if !defined(lrint)
#	define lrint (int64_t)rint
#endif
#endif
#ifdef NO_MEX
/* For testing, no data are actually touched and Matlab/Octave is not linked */
int ID_lhs = 100000;	/* We start output IDs at 100000 for fun */
int ID_rhs = 0;
#define mxArray void
#define mexErrMsgTxt(txt) fprintf (stderr, txt)
#endif

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

#define N_MEX_FIELDNAMES_GRID	22

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
	char    *fieldnames[N_MEX_FIELDNAMES_GRID];	/* this array contains the names of the fields of the output grid structure. */

	if (!G->data)
		mexErrMsgTxt ("GMTMEX_Get_Grid: programming error, output matrix G is empty\n");

	memset (fieldnames, 0, N_MEX_FIELDNAMES_GRID*sizeof (char *));
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
	grid_struct = mxCreateStructMatrix (1, 1, N_MEX_FIELDNAMES_GRID, (const char **)fieldnames );

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
	GMT_Report (API, GMT_MSG_NORMAL, "GMTMEX_Get_Text: INTERNAL ERROR: Handling of TEXTSET not implemented yet\n");
	return NULL;
}

#define N_MEX_FIELDNAMES_CPT	3

void * GMTMEX_Get_CPT (void *API, struct GMT_PALETTE *C)
{	/* Hook this Matlab CPT into the k'th output item */

	unsigned int k, j, n_colors;
	double *color = NULL, *alpha = NULL, *range = NULL;
	mxArray *mxcolormap = NULL, *mxalpha = NULL, *mxrange = NULL;
	mxArray *CPT_struct = NULL;
	char    *fieldnames[N_MEX_FIELDNAMES_CPT];	/* this array contains the names of the fields of the output grid structure. */

	if (!C->range)
		mexErrMsgTxt ("GMTMEX_Get_CPT: programming error, output CPT C is empty\n");

	memset (fieldnames, 0, N_MEX_FIELDNAMES_CPT*sizeof (char *));
	/* Return CPT via colormap, range, and alpha arrays in a struct */
	/* Create a Matlab struct for this CPT */
	fieldnames[0]  = mxstrdup ("colormap");
	fieldnames[1]  = mxstrdup ("range");
	fieldnames[2]  = mxstrdup ("alpha");
	CPT_struct = mxCreateStructMatrix (1, 1, N_MEX_FIELDNAMES_CPT, (const char **)fieldnames );

	n_colors = (C->is_continuous) ? C->n_colors + 1 : C->n_colors;
	mxcolormap = mxCreateNumericMatrix (n_colors, 3, mxDOUBLE_CLASS, mxREAL);
	color = mxGetPr (mxcolormap);
	mxalpha = mxCreateNumericMatrix (n_colors, 1, mxDOUBLE_CLASS, mxREAL);
	alpha = mxGetPr (mxalpha);
	mxrange = mxCreateNumericMatrix (2, 1, mxDOUBLE_CLASS, mxREAL);
	range = mxGetPr (mxrange);
	for (j = 0; j < C->n_colors; j++) {	/* Copy r/g/b from palette to Matlab array */
		for (k = 0; k < 3; k++) color[j+k*n_colors] = C->range[j].rgb_low[k];
		alpha[j] = C->range[j].rgb_low[3];
	}
	if (C->is_continuous) {	/* Add last color */
		for (k = 0; k < 3; k++) color[j+k*n_colors] = C->range[C->n_colors-1].rgb_high[k];
	}
	range[0] = C->range[0].z_low;
	range[1] = C->range[C->n_colors-1].z_high;
	
	mxSetField (CPT_struct, 0, "colormap", mxcolormap);
	mxSetField (CPT_struct, 0, "range", mxrange);
	mxSetField (CPT_struct, 0, "alpha", mxalpha);
	return (CPT_struct);
}

void * GMTMEX_Get_Image (void *API, struct GMT_MATRIX *M) {
	GMT_Report (API, GMT_MSG_NORMAL, "GMTMEX_Get_Image: INTERNAL ERROR: Handling of IMAGE not implemented yet\n");
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
			mexErrMsgTxt ("GMTMEX_grid_init: The input that was supposed to contain the Grid, is empty\n");
		if (!mxIsStruct (ptr))
			mexErrMsgTxt ("GMTMEX_grid_init: Expected a Grid structure for input\n");
		mx_ptr = mxGetField (ptr, 0, "inc");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_grid_init: Could not find inc array with Grid increments\n");
		inc = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "range");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_grid_init: Could not find range array for Grid range\n");
		range = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "registration");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_grid_init: Could not find registration array for Grid registration\n");
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
			mexErrMsgTxt ("GMTMEX_grid_init: Could not find data array for Grid\n");

		f = mxGetData (mx_ptr);
		for (gmt_ij = row = 0; row < G->header->ny; row++)
			for (col = 0; col < G->header->nx; col++, gmt_ij++)
				G->data [gmt_ij] = f[MEXG_IJ(G,row,col)];
		GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_grid_init: Allocated GMT Grid %lx\n", (long)G);
		GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_grid_init: Registered GMT Grid array %lx via memory reference from Matlab\n", (long)G->data);
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
		if (!mxIsNumeric (ptr)) mexErrMsgTxt ("GMTMEX_matrix_init: Expected a Matrix for input\n");
		dim[DIM_ROW] = mxGetM (ptr);	/* Number of rows */
		dim[DIM_COL] = mxGetN (ptr);	/* Number of columns */
		this_dim = dim;
	}
	/* Else there are no dimensions yet, this_dim is NULL and we are getting an empty container for output */
	if ((M = GMT_Create_Data (API, GMT_IS_MATRIX, GMT_IS_PLP, 0, this_dim, NULL, NULL, 0, 0, NULL)) == NULL)
		mexErrMsgTxt ("GMTMEX_matrix_init: Failure to alloc GMT source matrix\n");

	GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_matrix_init: Allocated GMT Matrix %lx\n", (long)M);
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
			mexErrMsgTxt ("GMTMEX_matrix_init: Unsupported Matlab data type in GMT matrix input.");
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

struct GMT_PALETTE *GMTMEX_CPT_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty CPT container to hold a GMT CPT.
 	 * If direction is GMT_IN then we are given a Matlab CPT and can determine its size, etc.
	 * If direction is GMT_OUT then we allocate an empty GMT CPT as a destination. */
	struct GMT_PALETTE *P = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		unsigned int k, j;
		uint64_t dim[1];
		mxArray *mx_ptr = NULL;
		double dz, *colormap = NULL, *range = NULL, *alpha = NULL;

		if (mxIsEmpty (ptr))
			mexErrMsgTxt ("GMTMEX_CPT_init: The input that was supposed to contain the CPT, is empty\n");
		if (!mxIsStruct (ptr))
			mexErrMsgTxt ("GMTMEX_CPT_init: Expected a CPT structure for input\n");
		mx_ptr = mxGetField (ptr, 0, "colormap");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_CPT_init: Could not find colormap array with CPT values\n");
		dim[0] = mxGetM (mx_ptr);	/* Number of rows */
		colormap = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "range");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_CPT_init: Could not find range array for CPT range\n");
		range = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "alpha");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_CPT_init: Could not find alpha array for CPT transparency\n");
		alpha = mxGetData (mx_ptr);
		if ((P = GMT_Create_Data (API, GMT_IS_CPT, GMT_IS_NONE, 0,
		                          dim, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_CPT_init: Failure to alloc GMT source CPT for input\n");
		dz = (range[1] - range[0]) / (P->n_colors + 1);
		for (j = 0; j < P->n_colors; j++) {
			for (k = 0; k < 3; k++)
				P->range[j].rgb_low[k] = colormap[j+k*dim[0]];
			P->range[j].rgb_low[3] = alpha[j];
			P->range[j].z_low = j * dz;
			P->range[j].z_high = (j+1) * dz;
		}
		GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_CPT_init: Allocated GMT CPT %lx\n", (long)P);
	}
	else {	/* Just allocate an empty container to hold an output grid (signal this by passing NULLs) */
		if ((P = GMT_Create_Data (API, GMT_IS_CPT, GMT_IS_NONE, 0,
                        NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_CPT_init: Failure to alloc GMT blank CPT container for holding output CPT\n");
	}
	return (P);
}

void * GMTMEX_Register_IO (void *API, unsigned int family, unsigned int geometry, unsigned int direction, const mxArray *ptr, int *ID)
{	/* Create the grid or matrix container, register it, and return the ID */
	void *obj = NULL;		/* Pointer to the container we created */
	*ID = GMT_NOTSET;

	switch (family) {
		case GMT_IS_GRID:
			/* Get an empty grid, and if input we and associate it with the Matlab grid pointer */
			obj = GMTMEX_grid_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_GRID, direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got Grid with Object ID %d\n", *ID);
			break;
		case GMT_IS_DATASET:
			/* Get a matrix container, and if input and associate it with the Matlab pointer */
			obj = GMTMEX_matrix_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_DATASET, direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got Matrix with Object ID %d\n", *ID);
			break;
		case GMT_IS_CPT:
			/* Get a CPT container, and if input and associate it with the Matlab CPT pointer */
			obj = GMTMEX_CPT_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_CPT, direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got CPT with Object ID %d\n", *ID);
			break;
		default:
			GMT_Report (API, GMT_MSG_NORMAL, "GMTMEX_Register_IO: Bad data type (%d)\n", family);
			break;
	}
	return (obj);
}
#endif
