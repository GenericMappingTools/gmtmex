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
 *      Contact info: www.soest.hawaii.edu/gmt
 *--------------------------------------------------------------------*/
/* GMT convenience functions used by MATLAB/OCTAVE mex functions.
 * All code that requires knowledge about MATLAB/OCTAVE functions is
 * found here.
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

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min macro */
#endif

enum MEX_dim {
	DIM_COL	= 0,	/* Holds the number of columns for vectors and x-nodes for matrix */
	DIM_ROW = 1};	/* Holds the number of rows for vectors and y-nodes for matrix */

/* Note: Wherever we say "Matlab" below we mean "Matlab of Octave" */

/* For the Mex interface we will wish to pass either filenames or mex variables via GMT command options.
 * We select a Matlab variable by suppying $ as the file name.  The parser will then find these $
 * arguments and replace them with references to mex variables via the GMT API mechanisms.
 * This requires us to know which options in a module may accept a file name.  As an example,
 * consider surface whose -Lu|l option may take a grid.  To pass a Matlab grid already in memory
 * we would use -Lu$ and give the grid as an argument to the module, e.g.,
 *    Z = gmt ('surface -R0/50/0/50 -I1 -V xyzfile -Lu$', uppergrid);
 * For each option that may take a file we need to know what kind of file and if this is input or output.
 * We encode this in a 3-character word XYZ, explained below.  Note that each module may
 * need several comma-separated XYZ words and these are returned as one string via GMT_Get_Moduleinfo.
 * The origin of these words are given by the THIS_MODULE_KEY in every module source code.
 *
 * X stands for the specific program option (e.g., L for -L, F for -F) or <,>
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
 *
 * All 5 GMT Resources are supported in this API, according to these rules:
 *  GMT_GRID:	Handled with a Matlab grid structure that holds, we use GMT's native GMT_GRID for the passing.
 *		  + Basic header array of length 9 [xmin, xmax, ymin, ymax, zmin, zmax, reg, xinc, yinc]
 *		  + The 2-D grid array (single precision)
 *		  + An x-array of coordinates
 *		  + An y-array of coordinates
 * GMT_DATASET: Handled with a Matlab matrix and we use GMT's native GMT_MATRIX for the passing.
 *		  + A 2-D matrix with rows and columns (double precision)
 * GMT_TEXTSET: Handled with a Matlab cell array and we use GMT's native GMT_TEXTSET for the passing.
 *		  + A 1-D cell array with one text record per cell
 * GMT_PALETTE: Handled with a Matlab structure and we use GMT's native GMT_PALETTE for the passing.
 *		  + colormap is the N*3 matrix for Matlab colormaps
 *		  + range is a 2-element array with zmin and zmax
 *		  + alpha is a N-element array with transparencies
 */

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
{	/* Given a GMT grid G, build a Matlab structure and assign the output components */

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

void * GMTMEX_Get_Dataset (void *API, struct GMT_MATRIX *M)
{	/* Given a GMT dataset via a matrix, build a Matlab matrix and assign values */
	unsigned int row, col;
	uint64_t gmt_ij, mex_ij;
	/* Create a 2-D Matlab double matrix of correct size */
	mxArray *P = mxCreateNumericMatrix (M->n_rows, M->n_columns, mxDOUBLE_CLASS, mxREAL);
	double *d = mxGetData (P);
	/* Duplicate the double GMT data matrix into the matlab double array */
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

void * GMTMEX_Get_Textset (void *API, struct GMT_TEXTSET *T)
{	/* Given a GMT textset T, build a Matlab cell array and assign values */
	uint64_t seg, row, k;
	mxArray *C = NULL, *p = NULL;
	struct GMT_TEXTSEGMENT *S = NULL;
	char text[BUFSIZ] = {""};
	
	if (T == NULL || !T->table)
		mexErrMsgTxt ("GMTMEX_Get_Textset: programming error, output textset T is NULL or empty\n");
	/* Create a cell array to hold all records */
	k = T->n_records;
	if (T->table[0]->n_segments > 1) k += T->table[0]->n_segments;
	C = mxCreateCellMatrix (k, 1);
	/* There is only one table when used in the external API, but it may have many segments.
	 * The segment information is lost when returned to Matlab */
	for (seg = k = 0; seg < T->table[0]->n_segments; seg++) {
		S = T->table[0]->segment[seg];
		if (T->table[0]->n_segments > 1) {
			sprintf (text, "> %s", S->header);
			p = mxCreateString (text);
			mxSetCell (C, k++, p);
		}
		for (row = 0; row <S->n_rows; row++, k++) {
			p = mxCreateString (S->record[row]);
			mxSetCell (C, k, p);
		}
	}
	return C;
}

void GMTMEX_Free_Textset (void *API, struct GMT_TEXTSET *T)
{
	/* Because of Windows DLL Hell we have to free those strdup'ed strings
	 * done in GMTMEX_init_text here instead of in GMT_Destroy_Data.
	 */
	uint64_t seg, row, k;
	struct GMT_TEXTSEGMENT *S = NULL;

	if (T == NULL || !T->table)
		mexErrMsgTxt ("GMTMEX_Get_Textset: programming error, textset T is NULL or empty\n");
	for (seg = k = 0; seg < T->table[0]->n_segments; seg++) {
		S = T->table[0]->segment[seg];
		for (row = 0; row <S->n_rows; row++) {
			free (S->record[row]);
			S->record[row] = NULL;
		}
	}
}

#define N_MEX_FIELDNAMES_CPT	3

void * GMTMEX_Get_CPT (void *API, struct GMT_PALETTE *C)
{	/* Given a GMT CPT C, build a Matlab structure and assign values */

	unsigned int k, j, n_colors;
	double *color = NULL, *alpha = NULL, *range = NULL;
	mxArray *mxcolormap = NULL, *mxalpha = NULL, *mxrange = NULL;
	mxArray *CPT_struct = NULL;
	char *fieldnames[N_MEX_FIELDNAMES_CPT];	/* Array with the names of the fields of the output grid structure. */

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
		alpha[j] = C->range[j].rgb_low[3];
	}
	range[0] = C->range[0].z_low;
	range[1] = C->range[C->n_colors-1].z_high;
	
	mxSetField (CPT_struct, 0, "colormap", mxcolormap);
	mxSetField (CPT_struct, 0, "range", mxrange);
	mxSetField (CPT_struct, 0, "alpha", mxalpha);
	return (CPT_struct);
}

#define N_MEX_FIELDNAMES_IMAGE	24

void *GMTMEX_Get_Image (void *API, struct GMT_IMAGE *I) {
	int item, n;
	mwSize dim[3];
	unsigned int row, col;
	uint64_t gmt_ij, mex_ij;
	uint8_t *u = NULL, *alpha = NULL;
	double *d = NULL, *dptr = NULL, *I_x = NULL, *I_y = NULL, *x = NULL, *y = NULL;
	double *color = NULL;
	mxArray *mxImg = NULL, *mx_x = NULL, *mx_y= NULL, *mxalpha = NULL, *mxcolormap = NULL;
	mxArray *mxProjectionRef = NULL;
	mxArray *mxHeader = NULL, *mxtmp = NULL;
	mxArray *image_struct = NULL;
	char    *fieldnames[N_MEX_FIELDNAMES_IMAGE];	/* this array contains the names of the fields of the output grid structure. */

	if (!I->data)
		mexErrMsgTxt ("GMTMEX_Get_Image: programming error, output image I is empty\n");

	memset (fieldnames, 0, N_MEX_FIELDNAMES_IMAGE*sizeof (char *));
	/* Return umage via a uint8_t (mxUINT8_CLASS) matrix in a struct */
	/* Create a Matlab struct for this image */
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
	fieldnames[18] = mxstrdup ("image");
	fieldnames[19] = mxstrdup ("x_units");
	fieldnames[20] = mxstrdup ("y_units");
	fieldnames[21] = mxstrdup ("z_units");
	fieldnames[22] = mxstrdup ("colormap");
	fieldnames[23] = mxstrdup ("alpha");
	image_struct = mxCreateStructMatrix (1, 1, N_MEX_FIELDNAMES_IMAGE, (const char **)fieldnames );

	mxtmp = mxCreateString (I->header->ProjRefPROJ4);
	mxSetField (image_struct, 0, (const char *) "ProjectionRefPROJ4", mxtmp);

	mxtmp = mxCreateString (I->header->ProjRefWKT);
	mxSetField (image_struct, 0, (const char *) "ProjectionRefWKT", mxtmp);

	mxHeader = mxCreateNumericMatrix (1, 9, mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr (mxHeader);
	for (n = 0; n < 4; n++) dptr[n] = I->header->wesn[n];
	dptr[4] = I->header->z_min;	dptr[5] = I->header->z_max;
	dptr[6] = I->header->registration;
	for (n = 0; n < 2; n++) dptr[n+7] = I->header->inc[n];
	mxSetField (image_struct, 0, "hdr", mxHeader);

	dptr = mxGetPr(mxtmp = mxCreateNumericMatrix (1, 4, mxDOUBLE_CLASS, mxREAL));
	for (n = 0; n < 4; n++) dptr[n] = I->header->wesn[n];
	mxSetField (image_struct, 0, "range", mxtmp);

	dptr = mxGetPr(mxtmp = mxCreateNumericMatrix (1, 2, mxDOUBLE_CLASS, mxREAL));
	for (n = 0; n < 2; n++) dptr[n] = I->header->inc[n];
	mxSetField (image_struct, 0, "inc", mxtmp);

	mxtmp = mxCreateDoubleScalar ((double)I->header->ny);
	mxSetField (image_struct, 0, (const char *) "n_rows", mxtmp);

	mxtmp = mxCreateDoubleScalar ((double)I->header->nx);
	mxSetField (image_struct, 0, (const char *) "n_columns", mxtmp);

	dptr = mxGetPr(mxtmp = mxCreateNumericMatrix (1, 2, mxDOUBLE_CLASS, mxREAL));
	dptr[0] = I->header->z_min;	dptr[1] = I->header->z_max;
	mxSetField (image_struct, 0, "MinMax", mxtmp);

	mxtmp = mxCreateDoubleScalar ((double)I->header->nan_value);
	mxSetField (image_struct, 0, (const char *) "NoDataValue", mxtmp);

	mxtmp = mxCreateDoubleScalar ((double)I->header->n_bands);
	mxSetField (image_struct, 0, (const char *) "LayerCount", mxtmp);
	
	dptr = mxGetPr(mxtmp = mxCreateNumericMatrix (1, 2, mxDOUBLE_CLASS, mxREAL));
	dptr[0] = I->header->ny;	dptr[1] = I->header->nx;
	mxSetField (image_struct, 0, "dim", mxtmp);

	mxtmp = mxCreateDoubleScalar ((double)I->header->registration);
	mxSetField (image_struct, 0, (const char *) "registration", mxtmp);

	mxtmp = mxCreateString (I->header->title);
	mxSetField (image_struct, 0, (const char *) "title", mxtmp);

	mxtmp = mxCreateString (I->header->command);
	mxSetField (image_struct, 0, (const char *) "command", mxtmp);

	mxtmp = mxCreateString (I->header->remark);
	mxSetField (image_struct, 0, (const char *) "remark", mxtmp);

	mxtmp = mxCreateString ("uint8");
	mxSetField (image_struct, 0, (const char *) "DataType", mxtmp);

	mxtmp = mxCreateDoubleScalar ((double)MIN(3,I->header->n_bands));
	mxSetField (image_struct, 0, (const char *) "LayerCount", mxtmp);

	mxtmp = mxCreateString (I->header->x_units);
	mxSetField (image_struct, 0, (const char *) "x_units", mxtmp);

	mxtmp = mxCreateString (I->header->y_units);
	mxSetField (image_struct, 0, (const char *) "y_units", mxtmp);

	mxtmp = mxCreateString (I->header->z_units);
	mxSetField (image_struct, 0, (const char *) "z_units", mxtmp);

	if (I->ColorMap != NULL) {	/* Indexed image has a color map */
		mxcolormap = mxCreateNumericMatrix (256, 3, mxDOUBLE_CLASS, mxREAL);
		mxImg = mxCreateNumericMatrix (I->header->ny, I->header->nx, mxUINT8_CLASS, mxREAL);
		color = mxGetPr (mxcolormap);
		u = mxGetData (mxImg);
		for (n = 0; n < 4 * 256 && I->ColorMap[n] >= 0; n++) color[n] = (uint8_t)I->ColorMap[n];
		n /= 4;
		memcpy (u, I->data, I->header->nm * sizeof (uint8_t));
		mxSetField (image_struct, 0, "colormap", mxcolormap);
	}	
	else if (I->header->n_bands == 1) { /* gray image */
		mxImg = mxCreateNumericMatrix (I->header->ny, I->header->nx, mxUINT8_CLASS, mxREAL);
		u = mxGetData (mxImg);
		memcpy (u, I->data, I->header->nm * sizeof (uint8_t));
	}
	else if (I->header->n_bands == 3) { /* RGB image */
		mxImg = mxCreateNumericMatrix (I->header->ny, I->header->nx, mxUINT8_CLASS, mxREAL);
		u = mxGetData (mxImg);
		memcpy (u, I->data, 3 * I->header->nm * sizeof (uint8_t));
	}
	else if (I->header->n_bands == 4) { /* RGBA image, with a color map */
		dim[0] = I->header->ny;	dim[1] = I->header->nx; dim[2] = 3;
		mxImg = mxCreateNumericArray (3, dim, mxUINT8_CLASS, mxREAL);
		u = mxGetData (mxImg);
		mxalpha = mxCreateNumericMatrix (I->header->ny, I->header->nx, mxUINT8_CLASS, mxREAL);
		alpha = mxGetData (mxalpha);
		memcpy(u, I->data, 3 * I->header->nm * sizeof (uint8_t)); 
		memcpy(alpha, &(I->data)[3 * I->header->nm], I->header->nm * sizeof (uint8_t)); 
		/*
		for (n = 0; n < I->header->nm; n++) {
			memcpy (&u[3*n], &(I->data)[4*n], 3 * sizeof (uint8_t));
			alpha[n] = (uint8_t)I->data[4*n+3];
		}
		*/
		mxSetField (image_struct, 0, "alpha", mxalpha);
	}
	mxSetField (image_struct, 0, "image", mxImg);

	/* Also return the convenient x and y arrays */
	I_x = GMT_Get_Coord (API, GMT_IS_IMAGE, GMT_X, I);	/* Get array of x coordinates */
	I_y = GMT_Get_Coord (API, GMT_IS_IMAGE, GMT_Y, I);	/* Get array of y coordinates */
	mx_x = mxCreateNumericMatrix (1, I->header->nx, mxDOUBLE_CLASS, mxREAL);
	mx_y = mxCreateNumericMatrix (1, I->header->ny, mxDOUBLE_CLASS, mxREAL);
	x = mxGetData (mx_x);
	y = mxGetData (mx_y);
	memcpy (x, I_x, I->header->nx * sizeof (double));
	for (n = 0; n < I->header->ny; n++) y[I->header->ny-1-n] = I_y[n];	/* Must reverse the y-array */
	if (GMT_Destroy_Data (API, &I_x))
		mexPrintf("Warning: Failure to delete I_x (x coordinate vector)\n");
	if (GMT_Destroy_Data (API, &I_y))
		mexPrintf("Warning: Failure to delete I_y (y coordinate vector)\n");
	mxSetField (image_struct, 0, "x", mx_x);
	mxSetField (image_struct, 0, "y", mx_y);
	return (image_struct);
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

struct GMT_IMAGE *GMTMEX_image_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Image container to hold a GMT image.
 	 * If direction is GMT_IN then we are given a Matlab image and can determine its size, etc.
	 * If direction is GMT_OUT then we allocate an empty GMT image as a destination. */
	unsigned int row, col;
	uint64_t gmt_ij;
	struct GMT_IMAGE *I = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		mxArray *mx_ptr = NULL;
		double *inc = NULL, *range = NULL, *reg = NULL, *MinMax = NULL;
		float *f = NULL;

		if (mxIsEmpty (ptr))
			mexErrMsgTxt ("GMTMEX_image_init: The input that was supposed to contain the Image, is empty\n");
		if (!mxIsStruct (ptr))
			mexErrMsgTxt ("GMTMEX_image_init: Expected a Image structure for input\n");
		mx_ptr = mxGetField (ptr, 0, "inc");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_image_init: Could not find inc array with Image increments\n");
		inc = mxGetData (mx_ptr);
		mx_ptr = mxGetField (ptr, 0, "range");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_image_init: Could not find range array for Image range\n");
		range = mxGetData (mx_ptr);
		if ((I = GMT_Create_Data (API, GMT_IS_IMAGE, GMT_IS_SURFACE, GMT_GRID_ALL,
			NULL, range, inc, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_image_init: Failure to alloc GMT source image for input\n");

		mx_ptr = mxGetField (ptr, 0, "MinMax");
		if (mx_ptr != NULL) {	/* Because we sent a NULL instead of the data array, z_min, z_max are not known. Use those from ptr */
			MinMax = mxGetData (mx_ptr);
			I->header->z_min = MinMax[0];
			I->header->z_max = MinMax[1];
		}

		mx_ptr = mxGetField (ptr, 0, "z");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("GMTMEX_image_init: Could not find data array for Grid\n");

		f = mxGetData (mx_ptr);
		for (gmt_ij = row = 0; row < I->header->ny; row++)
			for (col = 0; col < I->header->nx; col++, gmt_ij++)
				I->data [gmt_ij] = f[MEXG_IJ(I,row,col)];
		GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_image_init: Allocated GMT Image %lx\n", (long)I);
		GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_image_init: Registered GMT Image array %lx via memory reference from Matlab\n", (long)I->data);
	}
	else {	/* Just allocate an empty container to hold an output grid (signal this by passing NULLs) */
		if ((I = GMT_Create_Data (API, GMT_IS_IMAGE, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY,
                        NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_image_init: Failure to alloc GMT blank image container for holding output image\n");
	}
	return (I);
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
		dim[0] = mxGetM (mx_ptr);	/* Number of rows minus one since continous CPT */
		if (dim[0] < 1)
			mexErrMsgTxt ("GMTMEX_CPT_init: Colormap array has no CPT values\n");
		dim[0]--;	/* The number of CPT slices is one less since colormap is continuous */
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
		dz = (range[1] - range[0]) / P->n_colors;
		for (j = 0; j < P->n_colors; j++) {	/* OK to access j+1'th elemenent since length of colormap is P->n_colors+1 */
			for (k = 0; k < 3; k++) {
				P->range[j].rgb_low[k]  = colormap[j+k*dim[0]];
				P->range[j].rgb_high[k] = colormap[(j+1)+k*dim[0]];
			}
			P->range[j].rgb_low[3] = alpha[j];
			P->range[j].rgb_high[3] = alpha[j+1];
			P->range[j].z_low = range[0] + j * dz;
			P->range[j].z_high = P->range[j].z_low + dz;
		}
		P->is_continuous = 1;
		GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_CPT_init: Allocated GMT CPT %lx\n", (long)P);
	}
	else {	/* Just allocate an empty container to hold an output grid (signal this by passing NULLs) */
		if ((P = GMT_Create_Data (API, GMT_IS_CPT, GMT_IS_NONE, 0,
                        NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_CPT_init: Failure to alloc GMT blank CPT container for holding output CPT\n");
	}
	return (P);
}

struct GMT_TEXTSET *GMTMEX_Text_init (void *API, unsigned int direction, const mxArray *ptr)
{	/* Used to Create an empty Textset container to hold a GMT TEXTSET.
 	 * If direction is GMT_IN then we are given a Matlab cell array and can determine its size, etc.
	 * If direction is GMT_OUT then we allocate an empty GMT TEXTSET as a destination. */
	struct GMT_TEXTSET *T = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		uint64_t rec, dim[3] = { 1, 1, 0};
		mxArray *mx_ptr = NULL;
		char *txt = NULL;
		struct GMT_TEXTSEGMENT *S = NULL;

		if (mxIsEmpty (ptr))
			mexErrMsgTxt ("GMTMEX_Text_init: The input that was supposed to contain the Cell array is empty\n");
		if (!mxIsCell (ptr))
			mexErrMsgTxt ("GMTMEX_Text_init: Expected a Cell array for input\n");
		dim[GMT_ROW] = mxGetM (ptr);	/* Number of records */
		if (dim[GMT_ROW] == 1) {	/* Check if we got a transpose arrangement or just one record */
			rec = mxGetN (ptr);	/* Also possibly number of records */
			if (rec > 1) dim[GMT_ROW] = rec;	/* User gave row-vector of cells */
		}
		if ((T = GMT_Create_Data (API, GMT_IS_TEXTSET, GMT_IS_NONE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_Text_init: Failure to alloc GMT source TEXTSET for input\n");
		S = T->table[0]->segment[0];	/* Only one segment coming from Matlab */
		S->n_rows = dim[GMT_ROW];
		for (rec = 0; rec < S->n_rows; rec++) {
			mx_ptr = mxGetCell (ptr, rec);
			txt = mxArrayToString (mx_ptr);
			S->record[rec] = strdup (txt);
		}
		T->n_records = T->table[0]->n_records = S->n_rows;
		GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Text_init: Allocated GMT TEXTSET %lx\n", (long)T);
	}
	else {	/* Just allocate an empty container to hold an output grid (signal this by passing NULLs) */
		if ((T = GMT_Create_Data (API, GMT_IS_TEXTSET, GMT_IS_NONE, 0,
                        NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("GMTMEX_Text_init: Failure to alloc GMT blank TEXTSET container for holding output TEXT\n");
	}
	return (T);
}

void * GMTMEX_Register_IO (void *API, unsigned int family, unsigned int geometry, unsigned int direction, const mxArray *ptr, int *ID)
{	/* Create the grid or matrix container, register it, and return the ID */
	void *obj = NULL;		/* Pointer to the container we created */
	*ID = GMT_NOTSET;

	switch (family) {
		case GMT_IS_GRID:
			/* Get an empty grid, and if input we associate it with the Matlab grid pointer */
			obj = GMTMEX_grid_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_GRID, direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got Grid with Object ID %d\n", *ID);
			break;
		case GMT_IS_IMAGE:
			/* Get an empty image, and if input we associate it with the Matlab image pointer */
			obj = GMTMEX_image_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_IMAGE, direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got Image with Object ID %d\n", *ID);
			break;
		case GMT_IS_DATASET:
			/* Get a matrix container, and if input we associate it with the Matlab pointer */
			obj = GMTMEX_matrix_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_DATASET, direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got Matrix with Object ID %d\n", *ID);
			break;
		case GMT_IS_TEXTSET:
			/* Get a TEXTSET container, and if input we associate it with the Matlab pointer */
			obj = GMTMEX_Text_init (API, direction, ptr);
			*ID = GMT_Get_ID (API, GMT_IS_TEXTSET, direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got TEXTSET with Object ID %d\n", *ID);
			break;
		case GMT_IS_CPT:
			/* Get a CPT container, and if input we associate it with the Matlab CPT pointer */
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
