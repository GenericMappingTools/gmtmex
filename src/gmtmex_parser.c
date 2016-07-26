/*
 *	$Id$
 *
 *	Copyright (c) 1991-2016 by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis and F. Wobbe
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

enum MEX_dim {
	DIM_COL	= 0,	/* Holds the number of columns for vectors and x-nodes for matrix */
	DIM_ROW = 1};	/* Holds the number of rows for vectors and y-nodes for matrix */

/* Note: Wherever we say "MATLAB" below we mean "MATLAB or Octave"
 *
 * Reading and writing occur inside gmt_api.c via the standard GMT containers.
 * The packing up of GMT structures into MATLAB structures and vice versa happens after getting the
 * results out of the GMT API and before passing them back to MATLAB.
 *
 * All 6 GMT Resources are supported in this API, according to these rules:
 *  GMT_GRID:	Handled with a MATLAB grid structure and we use GMT's GMT_GRID for the passing
 *		  + Basic header array of length 9 [xmin, xmax, ymin, ymax, zmin, zmax, reg, xinc, yinc]
 *		  + The 2-D grid array z (single precision)
 *		  + An x-array of coordinates
 *		  + An y-array of coordinates
 *		  + Various Proj4 strings
 * GMT_DATASET: Handled with an array of MATLAB data structures and we use GMT's GMT_DATASET for passing in and out of GMT.
 *		Each MEX structure contain data for one segment:
 *		  + A text string for the segment header
 *		  + A cell array with strings, [] since no text in datasets
 *		  + A 2-D matrix with rows and columns (double precision)
 *		  + First segment may also have dataset comment and proj4/wtk strings
 * GMT_TEXTSET: Handled with an array of MATLAB data structures and we use GMT's GMT_TEXTSET for passing in and out of GMT.
 *		Each MEX structure contain data for one segment:
 *		  + A text string for the segment header
 *		  + A 2-D matrix with rows and columns (double precision) if there is numerical data in the first columns (or empty)
 *		  + A cell array with strings from columns that could not be deciphered as data.
 *		  + First segment may also have textset comment and proj4/wtk strings
 * GMT_PALETTE: Handled with a MATLAB structure and we use GMT's native GMT_PALETTE for the passing.
 *		  + colormap is the N*3 matrix for MATLAB colormaps
 *		  + range is a N-element array with z-values at color changes
 *		  + alpha is a N-element array with transparencies
 *		  + minmax is a 2-element array with zmin and zmax
 *		  + bnf is a 3x3-element matrix with back,fore,NaN colors
 *		  + depth is a 1-element matrix with color depth (1, 8, 24 bits)
 *		  + hinge is a 1-element matrix with hinge z-value (or NaN)
 *		  + cpt is a N*6-element matrix with original CPT slice values
 *		  + comment holds any Palette comments
 * GMT_POSTSCRIPT: Handled with a MATLAB structure and we use GMT's native GMT_POSTSCRIPT for the passing.
 *		  + postscript is the single string with all the PostScript code
 *		  + length is the number of bytes in the string
 *		  + mode is the overlay/trailer indicator
 *		  + comment holds any PostScript comments
 */

static char *mxstrdup (const char *s) {
	/* A strdup replacement to be used in Mexs to avoid memory leaks since the MATLAB
	   memory management will take care to free the memory allocated by this function */
	char *d = mxMalloc (strlen (s) + 1);
	if (d == NULL) return NULL;
	strcpy (d, s);
	return d;
}

int GMTMEX_print_func (FILE *fp, const char *message) {
	/* Replacement for GMT's gmt_print_func.  It is being used indirectly via
	 * API->print_func.  Purpose of this is to allow MATLAB (which cannot use
	 * printf) to reset API->print_func to this function via GMT_Create_Session.
	 * This allows GMT's errors and warnings to appear in MATLAB. */

	mexPrintf (message);
	return 0;
}

static int gmtmex_getMNK (const mxArray *p, int which) {
	/* Get number of columns or number of bands of a mxArray.
	   which = 0 to inquire n_rows
	         = 1 to inquire n_columns
	         = 2 to inquire n_bands
	         = ? ERROR
	*/
	int nx, ny, nBands, nDims;
	const mwSize *dim_array = NULL;

	nDims     = mxGetNumberOfDimensions(p);
	dim_array = mxGetDimensions(p);
	ny = dim_array[0];
	nx = dim_array[1];
	nBands = dim_array[2];
	if (nDims == 2) 	/* Otherwise it would stay undefined */
		nBands = 1;

	if (which == 0)
		return ny;
	else if (which == 1)
		return nx;
	else if (which == 2)
		return nBands;
	else
		mexErrMsgTxt("gmtmex_getMNK: Bad dimension number!");
	return -1;
}

static void gmtmex_quit_if_missing (const char *function, const char *field) {
	char buffer[128] = {""};
	sprintf (buffer , "%s: Could not find structure field %s\n", function, field);
	mexErrMsgTxt (buffer);
}

void *GMTMEX_Get_Grid (void *API, struct GMT_GRID *G) {
	/* Given an incoming GMT grid G, build a MATLAB structure and assign the output components.
 	 * Note: Incoming GMT grid has standard padding while MATLAB grid has none. */

	unsigned int k;
	uint64_t row, col, gmt_ij;
	float  *f = NULL;
	double *d = NULL, *G_x = NULL, *G_y = NULL, *x = NULL, *y = NULL;
	mxArray *G_struct = NULL, *mxptr[N_MEX_FIELDNAMES_GRID];

	if (!G->data)	/* Safety valve */
		mexErrMsgTxt ("GMTMEX_Get_Grid: programming error, output matrix G is empty\n");

	/* Create a MATLAB struct to hold this grid [matrix will be a float (mxSINGLE_CLASS)]. */
	G_struct = mxCreateStructMatrix (1, 1, N_MEX_FIELDNAMES_GRID, GMTMEX_fieldname_grid);

	/* Get pointers and populate structure from the information in G */
	mxptr[0]  = mxCreateNumericMatrix (G->header->n_rows, G->header->n_columns, mxSINGLE_CLASS, mxREAL);
	mxptr[1]  = mxCreateNumericMatrix (1, G->header->n_columns, mxDOUBLE_CLASS, mxREAL);
	mxptr[2]  = mxCreateNumericMatrix (1, G->header->n_rows,    mxDOUBLE_CLASS, mxREAL);
	mxptr[3]  = mxCreateNumericMatrix (1, 6, mxDOUBLE_CLASS, mxREAL);
	mxptr[4]  = mxCreateNumericMatrix (1, 2, mxDOUBLE_CLASS, mxREAL);
	mxptr[5]  = mxCreateDoubleScalar ((double)G->header->registration);
	mxptr[6]  = mxCreateDoubleScalar ((double)G->header->nan_value);
	mxptr[7]  = mxCreateString (G->header->title);
	mxptr[8]  = mxCreateString (G->header->remark);
	mxptr[9]  = mxCreateString (G->header->command);
	mxptr[10] = mxCreateString ("float32");
	mxptr[11] = mxCreateString (G->header->x_units);
	mxptr[12] = mxCreateString (G->header->y_units);
	mxptr[13] = mxCreateString (G->header->z_units);
	mxptr[14] = mxCreateString (G->header->ProjRefPROJ4);
	mxptr[15] = mxCreateString (G->header->ProjRefWKT);

	d = mxGetPr (mxptr[3]);	/* Range */
	for (k = 0; k < 4; k++) d[k] = G->header->wesn[k];
	d[4] = G->header->z_min;	d[5] = G->header->z_max;

	d = mxGetPr (mxptr[4]);	/* Increments */
	for (k = 0; k < 2; k++) d[k] = G->header->inc[k];

	/* Load the real grd array into a float MATLAB array by transposing
           from padded GMT grd format to unpadded MATLAB format */
	f = mxGetData (mxptr[0]);
	for (row = 0; row < G->header->n_rows; row++) {
		for (col = 0; col < G->header->n_columns; col++) {
			gmt_ij = GMT_IJP (G->header, row, col);
			f[MEXG_IJ(G,row,col)] = G->data[gmt_ij];
		}
	}

	/* Also return the convenient x and y arrays */
	G_x = GMT_Get_Coord (API, GMT_IS_GRID, GMT_X, G);	/* Get array of x coordinates */
	G_y = GMT_Get_Coord (API, GMT_IS_GRID, GMT_Y, G);	/* Get array of y coordinates */
	x = mxGetData (mxptr[1]);
	y = mxGetData (mxptr[2]);
	memcpy (x, G_x, G->header->n_columns * sizeof (double));
	for (k = 0; k < G->header->n_rows; k++)
		y[G->header->n_rows-1-k] = G_y[k];	/* Must reverse the y-array */
	if (GMT_Destroy_Data (API, &G_x))
		mexPrintf("Warning: Failure to delete G_x (x coordinate vector)\n");
	if (GMT_Destroy_Data (API, &G_y))
		mexPrintf("Warning: Failure to delete G_y (y coordinate vector)\n");
	for (k = 0; k < N_MEX_FIELDNAMES_GRID; k++)
		mxSetField (G_struct, 0, GMTMEX_fieldname_grid[k], mxptr[k]);
	return (G_struct);
}

void *GMTMEX_Get_Dataset (void *API, struct GMT_DATASET *D) {
	/* Given a GMT DATASET D, build a MATLAB array of segment structure and assign values.
	 * Each segment will have 6 items:
	 * header:	Text string with the segment header (could be empty)
	 * data:	Matrix with the data for this segment (n_rows by n_columns)
	 * text:	Empty cell array (since datasets have no text)
	 * comment:	Cell array with any comments
	 * proj4:	String with any proj4 information
	 * wkt:		String with any WKT information
	 */

	int n_headers;
	uint64_t tbl, seg, seg_out, col, start, k;
	double *data = NULL;
	struct GMT_DATASEGMENT *S = NULL;
	mxArray *D_struct = NULL, *mxheader = NULL, *mxdata = NULL, *mxtext = NULL, *mxstring = NULL;

	if (D == NULL || D->n_records == 0)	/* Safety valve */
		mexErrMsgTxt ("GMTMEX_Get_Dataset: programming error, output DATASET D is empty\n");
	
	for (tbl = seg_out = 0; tbl < D->n_tables; tbl++)	/* Count non-zero segments */
		for (seg = 0; seg < D->table[tbl]->n_segments; seg++)
			if (D->table[tbl]->segment[seg]->n_rows)
				seg_out++;
	
	D_struct = mxCreateStructMatrix ((mwSize)seg_out, 1, N_MEX_FIELDNAMES_DATASET, GMTMEX_fieldname_dataset);

	n_headers = D->table[0]->n_headers;
	for (tbl = seg_out = 0; tbl < D->n_tables; tbl++) {
		for (seg = 0; seg < D->table[tbl]->n_segments; seg++) {
			S = D->table[tbl]->segment[seg];	/* Shorthand */
			if (S->n_rows == 0) continue;		/* Skip empty segments */
			mxheader = mxCreateString (S->header);
			mxtext   = mxCreateCellMatrix (0, 0);	/* Empty */
			mxdata   = mxCreateNumericMatrix ((mwSize)S->n_rows, (mwSize)S->n_columns, mxDOUBLE_CLASS, mxREAL);
			data      = mxGetPr (mxdata);
			for (col = start = 0; col < S->n_columns; col++, start += S->n_rows) /* Copy the data columns */
				memcpy (&data[start], S->data[col], S->n_rows * sizeof (double));
			mxSetField (D_struct, (mwSize)seg_out, "data", mxdata);
			mxSetField (D_struct, (mwSize)seg_out, "text", mxtext);
			mxSetField (D_struct, (mwSize)seg_out, "header", mxheader);
			mxtext = mxCreateCellMatrix (n_headers, n_headers ? 1 : 0);
			for (k = 0; k < n_headers; k++) {
				mxstring = mxCreateString (D->table[0]->header[k]);
				mxSetCell (mxtext, (int)k, mxstring);
			}
			mxSetField (D_struct, (mwSize)seg_out, "comment", mxtext);
			n_headers = 0;	/* No other segment will have a non-empty comment cell array */
			seg_out++;
		}
	}
	return (D_struct);
}

int scan_to_start_of_text (char *text, uint64_t n_col) {
	/* Find the start of the n'th column where any optional text begins */
	uint64_t col = 0;
	int k = 0;
	if (n_col == 0) return 0;
	while (text[k] && strchr (" \t", text[k])) k++; /* Scan pass leading whitespace */
	while (text[k] && col < n_col) {
		while (text[k] && !strchr (" ,;\t", text[k])) k++;	/* Scan past this item until next "white-space" */
		col++;
		while (text[k] && strchr (" ,;\t", text[k])) k++;	/* Scan past consecutive "white-space" */
	}
	return k;
}

void *GMTMEX_Get_Textset (void *API, struct GMT_TEXTSET *T) {
	/* Given a GMT GMT_TEXTSET T, build a MATLAB array of segment structure and assign values.
	 * Each segment will have 6 items:
	 * header:	Text string with the segment header (could be empty)
	 * data:	Matrix with any converted data for this segment (n_rows by n_columns)
	 * text:	Cell array with the text items
	 * comment:	Cell array with any comments
	 * proj4:	String with any proj4 information
	 * wkt:		String with any WKT information
	 */

	int n_headers;
	uint64_t tbl, seg, seg_out, col, row, start, k, n_columns = 0;
	unsigned int flag[3] = {GMT_LAX_CONVERSION, 0, 0};	/* We will try to convert most of the text to data */
	double *data = NULL;
	struct GMT_DATASET *D = NULL;
	struct GMT_DATASEGMENT *SD = NULL;
	struct GMT_TEXTSEGMENT *ST = NULL;
	mxArray *D_struct = NULL, *mxheader = NULL, *mxdata = NULL, *mxtext = NULL, *mxstring = NULL;

	if (T == NULL || T->n_records == 0)	/* Safety valve */
		mexErrMsgTxt ("GMTMEX_Get_Textset: programming error, output GMT_TEXTSET T is empty\n");

	if ((D = GMT_Convert_Data (API, T, GMT_IS_TEXTSET, NULL, GMT_IS_DATASET, flag)) != NULL) {	/* Success */
		SD = D->table[0]->segment[0];	/* Shorthand, now determine number of non-NaN columns from first row */
		for (col = 0; col < D->n_columns; col++)
			if (!mxIsNaN (SD->data[col][0])) n_columns++;
	}

	for (tbl = seg_out = 0; tbl < T->n_tables; tbl++)	/* Count the non-empty segments */
		for (seg = 0; seg < T->table[tbl]->n_segments; seg++)
			if (T->table[tbl]->segment[seg]->n_rows)
				seg_out++;

	D_struct = mxCreateStructMatrix ((mwSize)seg_out, 1, N_MEX_FIELDNAMES_DATASET, GMTMEX_fieldname_dataset);
	n_headers = T->table[0]->n_headers;
	for (tbl = seg_out = 0; tbl < T->n_tables; tbl++) {
		for (seg = 0; seg < T->table[tbl]->n_segments; seg++) {
			ST = T->table[tbl]->segment[seg];	/* Shorthand to current text segment */
			if (ST->n_rows == 0) continue;		/* Skip empty segments */
			mxheader = mxCreateString (ST->header);
			if (D) {	/* We have numerial data to consider */
				SD     = D->table[tbl]->segment[seg];	/* Shorthand to the corresponding data segment */
				mxdata = mxCreateNumericMatrix ((mwSize)ST->n_rows, (mwSize)n_columns, mxDOUBLE_CLASS, mxREAL);
				data   = mxGetPr (mxdata);
				for (col = start = 0; col < n_columns; col++, start += SD->n_rows) /* Copy the data columns */
					memcpy (&data[start], SD->data[col], SD->n_rows * sizeof (double));
			}
			else	/* No data matrix, create an empty one */
				mxdata   = mxCreateNumericMatrix (0, 0, mxDOUBLE_CLASS, mxREAL);
			mxtext = mxCreateCellMatrix ((mwSize)ST->n_rows, 1);	/* Create cell array for text */
			for (row = 0; row < ST->n_rows; row++) {
				start = scan_to_start_of_text (ST->data[row], n_columns);
				mxstring = mxCreateString (&ST->data[row][start]);
				mxSetCell (mxtext, (int)row, mxstring);
			}
			mxSetField (D_struct, (mwSize)seg_out, "data",   mxdata);
			mxSetField (D_struct, (mwSize)seg_out, "text",   mxtext);
			mxSetField (D_struct, (mwSize)seg_out, "header", mxheader);
			mxtext = mxCreateCellMatrix (n_headers, n_headers ? 1 : 0);
			for (k = 0; k < n_headers; k++) {
				mxstring = mxCreateString (T->table[0]->header[k]);
				mxSetCell (mxtext, (int)k, mxstring);
			}
			mxSetField (D_struct, (mwSize)seg_out, "comment", mxtext);
			n_headers = 0;	/* No other segment will have a non-empty comment cell array */
			seg_out++;	/* Go to next output segment */
		}
	}
	if (D && GMT_Destroy_Data (API, &D))
		mexPrintf("Warning: Failure to delete intermediate D in GMTMEX_Get_Textset\n");
	return (D_struct);
}

#if GMT_MINOR_VERSION > 2

void *GMTMEX_Get_Postscript (void *API, struct GMT_POSTSCRIPT *P) {
	/* Given a GMT GMT_POSTSCRIPT P, build a MATLAB array of segment structure and assign values.
	 * Each segment will have 4 items:
	 * postscript:	Text string with the entire PostScript plot
	 * length:	Byte length of postscript
	 * mode:	1 has header, 2 has trailer, 3 is complete
	 * comment:	Cell array with any comments
	 */
	uint64_t k, *length = NULL;
	unsigned int *mode = NULL;
	mxArray *P_struct = NULL, *mxptr[N_MEX_FIELDNAMES_PS], *mxstring = NULL;
	
	if (P == NULL || !P->data)	/* Safety valve */
		mexErrMsgTxt ("GMTMEX_Get_Postscript: programming error, input POSTSCRIPT struct P is NULL or data string is empty\n");

	/* Return PS with postscript and length in a struct */
	P_struct = mxCreateStructMatrix (1, 1, N_MEX_FIELDNAMES_PS, GMTMEX_fieldname_ps);

	mxptr[0] = mxCreateString (P->data);
	mxptr[1] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
	mxptr[2] = mxCreateNumericMatrix (1, 1, mxUINT32_CLASS, mxREAL);
	mxptr[3] = mxCreateCellMatrix (P->n_headers, P->n_headers ? 1 : 0);
	length   = (uint64_t *)mxGetData(mxptr[1]);
	mode     = (uint32_t *)mxGetData(mxptr[2]);
	
	length[0] = (uint64_t)P->n_bytes;	/* Set length of the PS string */
	mode[0]   = (uint32_t)P->mode;		/* Set mode of the PS string */
	
	for (k = 0; k < P->n_headers; k++) {
		mxstring = mxCreateString (P->header[k]);
		mxSetCell (mxptr[3], (int)k, mxstring);
	}
	
	for (k = 0; k < N_MEX_FIELDNAMES_PS; k++)
		mxSetField (P_struct, 0, GMTMEX_fieldname_ps[k], mxptr[k]);

	return P_struct;
}
#endif

void *GMTMEX_Get_CPT (void *API, struct GMT_PALETTE *C) {
	/* Given a GMT GMT_PALETTE C, build a MATLAB structure and assign values.
	 * Each segment will have 10 items:
	 * colormap:	Nx3 array of colors usable in Matlab' colormap
	 * alpha:	Nx1 array with transparency values
	 * range:	Nx1 arran with z-values at color changes
	 * minmax:	2x1 array with min/max zvalues
	 * bfn:		3x3 array with colors for background, forground, nan
	 * depth	Color depth 24, 8, 1
	 * hinge:	Z-value at discontinuous color break, or NaN
	 * cpt:		Nx6 full GMT CPT array
	 * model:	String with color model rgb, hsv, or cmyk [rgb]
	 * comment:	Cell array with any comments
	 *
	 * Limitation: MATLAB's colormap format can either hold discrete
	 * or continuous colormaps, but not a mixture of these, which GMT
	 * can do.  Thus, mixed-mode GMT cpts being used in MATLAB or passed
	 * out from MATLAB cannot represent these changes accurately. */

	unsigned int k, j, n_colors, *depth = NULL;
	double *color = NULL, *cpt = NULL, *alpha = NULL, *minmax = NULL, *range = NULL, *hinge = NULL, *bfn = NULL;
	mxArray *C_struct = NULL, *mxptr[N_MEX_FIELDNAMES_CPT], *mxstring = NULL;

	if (C == NULL || !C->data)	/* Safety valve */
		mexErrMsgTxt ("GMTMEX_Get_CPT: programming error, output CPT C is empty\n");

	/* Return CPT via colormap, range, and alpha arrays in a struct */
	/* Create a MATLAB struct for this CPT */
	C_struct = mxCreateStructMatrix (1, 1, N_MEX_FIELDNAMES_CPT, GMTMEX_fieldname_cpt);

	n_colors = (C->is_continuous) ? C->n_colors + 1 : C->n_colors;
	mxptr[0] = mxCreateNumericMatrix (n_colors, 3, mxDOUBLE_CLASS, mxREAL);
	mxptr[1] = mxCreateNumericMatrix (n_colors, 1, mxDOUBLE_CLASS, mxREAL);
	mxptr[2] = mxCreateNumericMatrix (C->n_colors, 2, mxDOUBLE_CLASS, mxREAL);
	mxptr[3] = mxCreateNumericMatrix (2, 1, mxDOUBLE_CLASS, mxREAL);
	mxptr[4] = mxCreateNumericMatrix (3, 3, mxDOUBLE_CLASS, mxREAL);
	mxptr[5] = mxCreateNumericMatrix (1, 1, mxUINT32_CLASS, mxREAL);
	mxptr[6] = mxCreateNumericMatrix (1, 1, mxDOUBLE_CLASS, mxREAL);
	mxptr[7] = mxCreateNumericMatrix (C->n_colors, 6, mxDOUBLE_CLASS, mxREAL);
	mxptr[8] = NULL;	/* Set below */
	mxptr[9] = mxCreateCellMatrix (C->n_headers, C->n_headers ? 1 : 0);
	
	color    = mxGetPr (mxptr[0]);
	alpha    = mxGetPr (mxptr[1]);
	range    = mxGetPr (mxptr[2]);
	minmax   = mxGetPr (mxptr[3]);
	bfn      = mxGetPr (mxptr[4]);
	depth    = (uint32_t *)mxGetData (mxptr[5]);
	hinge    = mxGetPr (mxptr[6]);
	cpt      = mxGetPr (mxptr[7]);
	depth[0] = (C->is_bw) ? 1 : ((C->is_gray) ? 8 : 24);
	hinge[0] = (C->has_hinge) ? C->hinge : mxGetNaN ();
	for (j = 0; j < 3; j++)	/* Copy r/g/b from palette bfn to MATLAB array */
		for (k = 0; k < 3; k++) bfn[j+3*k] = C->bfn[j].rgb[k];
	for (j = 0; j < C->n_colors; j++) {	/* Copy r/g/b from palette to MATLAB colormap and cpt */
		for (k = 0; k < 3; k++) {
			color[j+k*n_colors] = cpt[j+k*C->n_colors] = C->data[j].rgb_low[k];
			cpt[j+(k+3)*C->n_colors] = C->data[j].rgb_high[k];
		}
		alpha[j] = C->data[j].rgb_low[3];
		range[j] = C->data[j].z_low;
		range[j+C->n_colors] = C->data[j].z_high;
	}
	if (C->is_continuous) {	/* Add last color/alpha to colormap */
		for (k = 0; k < 3; k++) color[j+k*n_colors] = C->data[C->n_colors-1].rgb_high[k];
		alpha[j] = C->data[j].rgb_low[3];
	}
	minmax[0] = C->data[0].z_low;	/* Set min/max limits */
	minmax[1] = C->data[C->n_colors-1].z_high;
	if (C->n_headers) {
		for (k = 0; k < C->n_headers; k++) {
			mxstring = mxCreateString (C->header[k]);
			mxSetCell (mxptr[9], (int)k, mxstring);
		}
	}
	if (C->model & GMT_HSV)
		mxptr[8] = mxCreateString ("hsv");
	else if (C->model & GMT_CMYK)
		mxptr[8] = mxCreateString ("cmyk");
	else
		mxptr[8] = mxCreateString ("rgb");

	for (k = 0; k < N_MEX_FIELDNAMES_CPT; k++)	/* Update all fields */
		mxSetField (C_struct, 0, GMTMEX_fieldname_cpt[k], mxptr[k]);
	return (C_struct);
}

void *GMTMEX_Get_Image (void *API, struct GMT_IMAGE *I) {
	unsigned int k, row, col;
	mwSize   dim[3], kk, m;
	uint8_t *u = NULL, *alpha = NULL;
	double  *d = NULL, *I_x = NULL, *I_y = NULL, *x = NULL, *y = NULL, *color = NULL;
	mxArray *I_struct = NULL, *mxptr[N_MEX_FIELDNAMES_IMAGE];

	if (I == NULL || !I->data)	/* Safety valve */
		mexErrMsgTxt ("GMTMEX_Get_Image: programming error, output image I is empty\n");

	/* Return image via a uint8_t (mxUINT8_CLASS) matrix in a struct */
	/* Create a MATLAB struct for this image */
	I_struct = mxCreateStructMatrix (1, 1, N_MEX_FIELDNAMES_IMAGE, GMTMEX_fieldname_image);
	/* Create the various fields with information from I */
	mxptr[0]  = NULL;	/* Set below */
	mxptr[1]  = mxCreateNumericMatrix (1, I->header->n_columns, mxDOUBLE_CLASS, mxREAL);
	mxptr[2]  = mxCreateNumericMatrix (1, I->header->n_rows, mxDOUBLE_CLASS, mxREAL);
	mxptr[3]  = mxCreateNumericMatrix (1, 6, mxDOUBLE_CLASS, mxREAL);
	mxptr[4]  = mxCreateNumericMatrix (1, 2, mxDOUBLE_CLASS, mxREAL);
	mxptr[5]  = mxCreateDoubleScalar ((double)I->header->registration);
	mxptr[6]  = mxCreateDoubleScalar ((double)I->header->nan_value);	
	mxptr[7]  = mxCreateString (I->header->title);
	mxptr[8]  = mxCreateString (I->header->remark);
	mxptr[9]  = mxCreateString (I->header->command);
	mxptr[10] = mxCreateString ("uint8");
	mxptr[11] = mxCreateString (I->header->x_units);
	mxptr[12] = mxCreateString (I->header->y_units);
	mxptr[13] = mxCreateString (I->header->z_units);
	mxptr[14] = mxptr[15] = NULL;	/* Set below */
	mxptr[16] = (I->header->mem_layout[0]) ? mxCreateString(I->header->mem_layout) : mxCreateString ("TCBa");
	mxptr[17] = mxCreateString (I->header->ProjRefPROJ4);
	mxptr[18] = mxCreateString (I->header->ProjRefWKT);

	/* Fill in values */
	d = mxGetPr (mxptr[3]);	/* Range */
	for (k = 0; k < 4; k++) d[k] = I->header->wesn[k];
	d[4] = I->header->z_min;	d[5] = I->header->z_max;

	d = mxGetPr(mxptr[4]);	/* Increments */
	for (k = 0; k < 2; k++) d[k] = I->header->inc[k];

	if (I->colormap != NULL) {	/* Indexed image has a color map */
		mxptr[14] = mxCreateNumericMatrix (I->n_indexed_colors, 3, mxDOUBLE_CLASS, mxREAL);
		mxptr[0]  = mxCreateNumericMatrix (I->header->n_rows, I->header->n_columns, mxUINT8_CLASS, mxREAL);
		u     = mxGetData (mxptr[0]);
		color = mxGetPr (mxptr[14]);
		for (k = 0; k < 4 * (unsigned int)I->n_indexed_colors && I->colormap[k] >= 0; k++)
			color[k] = (uint8_t)I->colormap[k];
		k /= 4;
		memcpy (u, I->data, I->header->nm * sizeof (uint8_t));
	}	
	else if (I->header->n_bands == 1) {	/* gray image */
		mxptr[0] = mxCreateNumericMatrix (I->header->n_rows, I->header->n_columns, mxUINT8_CLASS, mxREAL);
		u = mxGetData (mxptr[0]);
		memcpy (u, I->data, I->header->nm * sizeof (uint8_t));
	}
	else if (I->header->n_bands == 3) {	/* RGB image */
		dim[0] = I->header->n_rows;	dim[1] = I->header->n_columns; dim[2] = 3;
		mxptr[0] = mxCreateNumericArray (3, dim, mxUINT8_CLASS, mxREAL);
		u = mxGetData (mxptr[0]);
		if (!strncmp(I->header->mem_layout, "TCBa", 4))
			memcpy (u, I->data, 3 * I->header->nm * sizeof (uint8_t));
		else if (!strncmp(I->header->mem_layout, "TRPa", 4)) {
			kk = 0;
			for (row = 0; row < I->header->n_rows; row++)
				for (col = 0; col < I->header->n_columns; col++)
					for (m = 0; m < 3; m++)
						u[row + col*I->header->n_rows + m*I->header->nm] = (uint8_t)I->data[k++];
			mxptr[16] = mxCreateString ("TCBa");	/* Because we just converted to it above */
		}
		if (I->alpha) {
			mxptr[15] = mxCreateNumericMatrix (I->header->n_rows, I->header->n_columns, mxUINT8_CLASS, mxREAL);
			alpha = mxGetData (mxptr[15]);
			memcpy (alpha, I->alpha, I->header->nm * sizeof (uint8_t)); 
		}
	}
	else if (I->header->n_bands == 4) {	/* RGBA image, with a color map */
		dim[0] = I->header->n_rows;	dim[1] = I->header->n_columns; dim[2] = 3;
		mxptr[0] = mxCreateNumericArray (3, dim, mxUINT8_CLASS, mxREAL);
		u = mxGetData (mxptr[0]);
		mxptr[15] = mxCreateNumericMatrix (I->header->n_rows, I->header->n_columns, mxUINT8_CLASS, mxREAL);
		alpha = mxGetData (mxptr[15]);
		memcpy (u, I->data, 3 * I->header->nm * sizeof (uint8_t)); 
		memcpy (alpha, &(I->data)[3 * I->header->nm], I->header->nm * sizeof (uint8_t)); 
		/*
		for (k = 0; k < I->header->nm; k++) {
			for (m = 0; m < 3; m++)
				u[k+m*I->header->nm] = (uint8_t)I->data[4*k+m];
			alpha[k] = (uint8_t)I->data[4*k+3];
		}
		*/
	}

	/* Also return the convenient x and y arrays */
	I_x = GMT_Get_Coord (API, GMT_IS_IMAGE, GMT_X, I);	/* Get array of x coordinates */
	I_y = GMT_Get_Coord (API, GMT_IS_IMAGE, GMT_Y, I);	/* Get array of y coordinates */
	x = mxGetData (mxptr[1]);
	y = mxGetData (mxptr[2]);
	memcpy (x, I_x, I->header->n_columns * sizeof (double));
	for (k = 0; k < I->header->n_rows; k++)	/* Must reverse the y-array */
		y[I->header->n_rows-1-k] = I_y[k];
	if (GMT_Destroy_Data (API, &I_x))
		mexPrintf("Warning: Failure to delete I_x (x coordinate vector)\n");
	if (GMT_Destroy_Data (API, &I_y))
		mexPrintf("Warning: Failure to delete I_y (y coordinate vector)\n");
	for (k = 0; k < N_MEX_FIELDNAMES_IMAGE; k++) {	/* Update fields */
		if (mxptr[k]) mxSetField (I_struct, 0, GMTMEX_fieldname_image[k], mxptr[k]);
	}
	return (I_struct);
}

static struct GMT_GRID *gmtmex_grid_init (void *API, unsigned int direction, unsigned int module_input, const mxArray *ptr) {
	/* Used to Create an empty Grid container to hold a GMT grid.
 	 * If direction is GMT_IN then we are given a MATLAB grid and can determine its size, etc.
	 * If direction is GMT_OUT then we allocate an empty GMT grid as a destination. */
	unsigned int row, col;
	uint64_t gmt_ij;
	struct GMT_GRID *G = NULL;

	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		unsigned int registration, flag = (module_input) ? GMT_VIA_MODULE_INPUT : 0;
		mxArray *mx_ptr = NULL, *mxGrid = NULL, *mxHdr = NULL;

		if (mxIsEmpty (ptr))
			mexErrMsgTxt ("gmtmex_grid_init: The input that was supposed to contain the Grid, is empty\n");
		if (!mxIsStruct (ptr)) {
			if (!mxIsCell (ptr))
				mexErrMsgTxt ("gmtmex_grid_init: Expected a Grid structure or Cell array for input\n");
			else {		/* Test that we have a {MxN,1x9} cell array */
				if (mxGetM(ptr) != 2 && mxGetN(ptr) != 2)
					mexErrMsgTxt ("gmtmex_grid_init: Cell array must contain two elements\n");
				else {
					mxGrid = mxGetCell(ptr, 0);
					mxHdr  = mxGetCell(ptr, 1);
					if (mxGetM(mxGrid) < 2 || mxGetN(mxGrid) < 2)
						mexErrMsgTxt ("gmtmex_grid_init: First element of grid's cell array must contain a decent matrix\n");
					if (mxGetM(mxHdr) != 1 || mxGetN(mxHdr) != 9)
						mexErrMsgTxt ("gmtmex_grid_init: grid's cell array second element must contain a 1x9 vector\n");
					if (!mxIsSingle(mxGrid) && !mxIsDouble(mxGrid))
						mexErrMsgTxt ("gmtmex_grid_init: grid's cell matrix must be either single or double.\n");
				}
			}
		}

		if (mxIsStruct(ptr)) {	/* Passed a regular MEX Grid structure */
			double *inc = NULL, *range = NULL, *reg = NULL;
			char x_unit[GMT_GRID_VARNAME_LEN80] = { "" }, y_unit[GMT_GRID_VARNAME_LEN80] = { "" },
			     z_unit[GMT_GRID_VARNAME_LEN80] = { "" };
			mx_ptr = mxGetField (ptr, 0, "inc");
			if (mx_ptr == NULL)
				mexErrMsgTxt ("gmtmex_grid_init: Could not find inc array with Grid increments\n");
			inc = mxGetData (mx_ptr);

			mx_ptr = mxGetField (ptr, 0, "range");
			if (mx_ptr == NULL)
				mexErrMsgTxt ("gmtmex_grid_init: Could not find range array for Grid range\n");
			range = mxGetData (mx_ptr);

			mxGrid = mxGetField(ptr, 0, "z");
			if (mxGrid == NULL)
				mexErrMsgTxt ("gmtmex_grid_init: Could not find data array for Grid\n");
			if (!mxIsSingle(mxGrid) && !mxIsDouble(mxGrid))
				mexErrMsgTxt ("gmtmex_grid_init: data array must be either single or double.\n");

			mx_ptr = mxGetField (ptr, 0, "registration");
			if (mx_ptr == NULL)
				mexErrMsgTxt ("gmtmex_grid_init: Could not find registration array for Grid registration\n");
			reg = mxGetData (mx_ptr);
			registration = (unsigned int)lrint(reg[0]);
			if ((G = GMT_Create_Data (API, GMT_IS_GRID|flag, GMT_IS_SURFACE, GMT_GRID_ALL,
			                          NULL, range, inc, registration, GMT_NOTSET, NULL)) == NULL)
				mexErrMsgTxt ("gmtmex_grid_init: Failure to alloc GMT source matrix for input\n");

			G->header->z_min = range[4];
			G->header->z_max = range[5];

			G->header->registration = registration;

			mx_ptr = mxGetField (ptr, 0, "nodata");
			if (mx_ptr != NULL)
				G->header->nan_value = *(float *)mxGetData (mx_ptr);

			mx_ptr = mxGetField (ptr, 0, "proj4");
			if (mx_ptr != NULL && mxGetN(mx_ptr) > 6) {		/* A true proj4 string will have at least this lenght */
				char *str = malloc(mxGetN(mx_ptr) + 1);
				mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr));
				G->header->ProjRefPROJ4 = GMT_Duplicate_String (API, str);
				free (str);
			}
			mx_ptr = mxGetField (ptr, 0, "wkt");
			if (mx_ptr != NULL && mxGetN(mx_ptr) > 20) {	/* A true WTT string will have more thna this lenght */ 
				char *str = malloc(mxGetN(mx_ptr) + 1);
				mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr) + 1);
				G->header->ProjRefWKT = GMT_Duplicate_String (API, str);
				free (str);
			}
			mx_ptr = mxGetField (ptr, 0, "title");
			if (mx_ptr != NULL) {
				char *str = malloc(mxGetN(mx_ptr) + 1);
				mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr) + 1);
				strncpy(G->header->title, str, GMT_GRID_VARNAME_LEN80 - 1);
				free (str);
			}
			mx_ptr = mxGetField (ptr, 0, "command");
			if (mx_ptr != NULL) {
				char *str = malloc(mxGetN(mx_ptr) + 1);
				mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr) + 1);
				strncpy(G->header->command, str, GMT_GRID_COMMAND_LEN320 - 1);
				free (str);
			}
			mx_ptr = mxGetField (ptr, 0, "comment");
			if (mx_ptr != NULL) {
				char *str = malloc(mxGetN(mx_ptr)+2);
				mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr));
				strncpy(G->header->remark, str, GMT_GRID_REMARK_LEN160 - 1);
				free (str);
			}
			mx_ptr = mxGetField (ptr, 0, "x_unit");
			if (mx_ptr != NULL) {
				mxGetString(mx_ptr, x_unit, (mwSize)mxGetN(mx_ptr) + 1);
				strncpy(G->header->x_units, x_unit, GMT_GRID_VARNAME_LEN80 - 1);
			}
			mx_ptr = mxGetField (ptr, 0, "y_unit");
			if (mx_ptr != NULL) {
				mxGetString(mx_ptr, y_unit, (mwSize)mxGetN(mx_ptr) + 1);
				strncpy(G->header->y_units, y_unit, GMT_GRID_VARNAME_LEN80 - 1);
			}
			mx_ptr = mxGetField (ptr, 0, "z_unit");
			if (mx_ptr != NULL) {
				mxGetString(mx_ptr, z_unit, (mwSize)mxGetN(mx_ptr) + 1);
				strncpy(G->header->z_units, z_unit, GMT_GRID_VARNAME_LEN80 - 1);
			}
		}
		else {	/* Passed header and grid separately */
			double *h = mxGetData(mxHdr);
			registration = (unsigned int)lrint(h[6]);
			if ((G = GMT_Create_Data (API, GMT_IS_GRID|flag, GMT_IS_SURFACE, GMT_GRID_ALL,
			                          NULL, h, &h[7], registration, GMT_NOTSET, NULL)) == NULL)
				mexErrMsgTxt ("gmtmex_grid_init: Failure to alloc GMT source matrix for input\n");
			G->header->z_min = h[4];
			G->header->z_max = h[5];
		}

		if (mxIsSingle(mxGrid)) {
			float *f4 = mxGetData(mxGrid);
			if (f4 == NULL)
				mexErrMsgTxt("gmtmex_grid_init: Grid pointer is NULL where it absolutely could not be.");
			for (row = 0; row < G->header->n_rows; row++) {
				for (col = 0; col < G->header->n_columns; col++) {
					gmt_ij = GMT_IJP (G->header, row, col);
					G->data[gmt_ij] = f4[MEXG_IJ(G,row,col)];
				}
			}
		}
		else {
			double *f8 = mxGetData(mxGrid);
			if (f8 == NULL)
				mexErrMsgTxt("gmtmex_grid_init: Grid pointer is NULL where it absolutely could not be.");
			for (row = 0; row < G->header->n_rows; row++) {
				for (col = 0; col < G->header->n_columns; col++) {
					gmt_ij = GMT_IJP (G->header, row, col);
					G->data[gmt_ij] = (float)f8[MEXG_IJ(G,row,col)];
				}
			}
		}
		GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_grid_init: Allocated GMT Grid %lx\n", (long)G);
		GMT_Report (API, GMT_MSG_DEBUG,
		            "gmtmex_grid_init: Registered GMT Grid array %lx via memory reference from MATLAB\n",
		            (long)G->data);
	}
	else {	/* Just allocate an empty container to hold an output grid (signal this by passing 0s and NULLs) */
		if ((G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, 0,
		                          NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_grid_init: Failure to alloc GMT blank grid container for holding output grid\n");
	}
	return (G);
}

static struct GMT_IMAGE *gmtmex_image_init (void *API, unsigned int direction, unsigned int module_input, const mxArray *ptr) {
	/* Used to Create an empty Image container to hold a GMT image.
 	 * If direction is GMT_IN then we are given a MATLAB image and can determine its size, etc.
	 * If direction is GMT_OUT then we allocate an empty GMT image as a destination. */
	struct GMT_IMAGE *I = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		uint64_t dim[3];
		unsigned int flag = (module_input) ? GMT_VIA_MODULE_INPUT : 0;
		char x_unit[GMT_GRID_VARNAME_LEN80] = { "" }, y_unit[GMT_GRID_VARNAME_LEN80] = { "" },
		     z_unit[GMT_GRID_VARNAME_LEN80] = { "" }, layout[4];
		double  *reg = NULL, *inc = NULL, *range = NULL;
		mxArray *mx_ptr = NULL;

		if (mxIsEmpty (ptr))
			mexErrMsgTxt ("gmtmex_image_init: The input that was supposed to contain the Image, is empty\n");

		if (!mxIsStruct (ptr))
			mexErrMsgTxt ("gmtmex_image_init: Expected a Image structure for input\n");

		mx_ptr = mxGetField (ptr, 0, "range");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("gmtmex_image_init: Could not find range array for Image range\n");
		range = mxGetData (mx_ptr);

		mx_ptr = mxGetField (ptr, 0, "inc");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("gmtmex_image_init: Could not find inc array with Image increments\n");
		inc = mxGetData (mx_ptr);

		mx_ptr = mxGetField(ptr, 0, "registration");
		if (mx_ptr == NULL)
			mexErrMsgTxt("gmtmex_image_init: Could not find registration info in Image struct\n");
		reg = mxGetData(mx_ptr);

		mx_ptr = mxGetField (ptr, 0, "image");
		if (mx_ptr == NULL)
			mexErrMsgTxt ("gmtmex_image_init: Could not find data array for Image\n");

		if (!mxIsUint8(mx_ptr))
			mexErrMsgTxt("gmtmex_image_init: The only data type supported by now is UInt8, and this image is not.\n");

		dim[0] = gmtmex_getMNK (mx_ptr, 1);	dim[1] = gmtmex_getMNK (mx_ptr, 0);	dim[2] = gmtmex_getMNK (mx_ptr, 2);
		if ((I = GMT_Create_Data (API, GMT_IS_IMAGE|flag, GMT_IS_SURFACE, GMT_GRID_HEADER_ONLY, dim,
			                      range, inc, (unsigned int)reg[0], 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_image_init: Failure to alloc GMT source image for input\n");

		I->data = (unsigned char *)mxGetData (mx_ptr);				/* Send in the Matlab owned memory. */
		I->alloc_mode = GMT_ALLOC_EXTERNALLY;

/*
		memcpy (I->data, (unsigned char *)mxGetData (mx_ptr), I->header->nm * I->header->n_bands * sizeof (char));
		for (row = 0; row < I->header->n_rows; row++) {
			for (col = 0; col < I->header->n_columns; col++) {
				gmt_ij = GMT_IJP (I->header, row, col);
				I->data [gmt_ij] = f[MEXG_IJ(I,row,col)];
			}
		}
*/
		mx_ptr = mxGetField (ptr, 0, "alpha");
		I->alpha = NULL;
		if (mx_ptr != NULL) {
			if (mxGetNumberOfDimensions(mx_ptr) == 2)
				I->alpha = (unsigned char *)mxGetData (mx_ptr);		/* Send in the Matlab owned memory. */
		}

		I->header->z_min = range[4];
		I->header->z_max = range[5];

		mx_ptr = mxGetField(ptr, 0, "x");
		if (mx_ptr == NULL)
			mexErrMsgTxt("gmtmex_image_init: Could not find x-coords vector for Image\n");
		I->x = mxGetData(mx_ptr);

		mx_ptr = mxGetField(ptr, 0, "y");
		if (mx_ptr == NULL)
			mexErrMsgTxt("gmtmex_image_init: Could not find y-coords vector for Image\n");
		I->y = mxGetData(mx_ptr);

		mx_ptr = mxGetField (ptr, 0, "nodata");
		if (mx_ptr != NULL)
			I->header->nan_value = *(float *)mxGetData (mx_ptr);

		mx_ptr = mxGetField (ptr, 0, "proj4");
		if (mx_ptr != NULL && mxGetN(mx_ptr) > 6) {		/* A true proj4 string will have at least this lenght */
			char *str = malloc(mxGetN(mx_ptr) + 1);
			mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr) + 1);
			I->header->ProjRefPROJ4 = GMT_Duplicate_String (API, str);
			free (str);
		}
		mx_ptr = mxGetField (ptr, 0, "wkt");
		if (mx_ptr != NULL && mxGetN(mx_ptr) > 20) {	/* A true WTT string will have more thna this lenght */ 
			char *str = malloc(mxGetN(mx_ptr) + 1);
			mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr) + 1);
			I->header->ProjRefWKT = GMT_Duplicate_String (API, str);
			free (str);
		}

		mx_ptr = mxGetField (ptr, 0, "title");
		if (mx_ptr != NULL) {
			char *str = malloc(mxGetN(mx_ptr) + 1);
			mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr) + 1);
			strncpy(I->header->title, str, GMT_GRID_VARNAME_LEN80 - 1);
			free (str);
		}
		mx_ptr = mxGetField (ptr, 0, "command");
		if (mx_ptr != NULL) {
			char *str = malloc(mxGetN(mx_ptr) + 1);
			mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr) + 1);
			strncpy(I->header->command, str, GMT_GRID_COMMAND_LEN320 - 1);
			free (str);
		}
		mx_ptr = mxGetField (ptr, 0, "comment");
		if (mx_ptr != NULL) {
			char *str = malloc(mxGetN(mx_ptr) + 1);
			mxGetString(mx_ptr, str, (mwSize)mxGetN(mx_ptr) + 1);
			strncpy(I->header->remark, str, GMT_GRID_REMARK_LEN160 - 1);
			free (str);
		}
		mx_ptr = mxGetField (ptr, 0, "x_unit");
		if (mx_ptr != NULL) {
			mxGetString(mx_ptr, x_unit, (mwSize)mxGetN(mx_ptr));
			strncpy(I->header->x_units, x_unit, GMT_GRID_VARNAME_LEN80 - 1);
		}
		mx_ptr = mxGetField (ptr, 0, "y_unit");
		if (mx_ptr != NULL) {
			mxGetString(mx_ptr, y_unit, (mwSize)mxGetN(mx_ptr));
			strncpy(I->header->y_units, y_unit, GMT_GRID_VARNAME_LEN80 - 1);
		}
		mx_ptr = mxGetField (ptr, 0, "z_unit");
		if (mx_ptr != NULL) {
			mxGetString(mx_ptr, z_unit, (mwSize)mxGetN(mx_ptr));
			strncpy(I->header->z_units, z_unit, GMT_GRID_VARNAME_LEN80 - 1);
		}
		mx_ptr = mxGetField (ptr, 0, "layout");
		if (mx_ptr != NULL) {
			mxGetString(mx_ptr, layout, (mwSize)mxGetN(mx_ptr));
			strncpy(I->header->mem_layout, layout, 4);
		}
		else
			strncpy(I->header->mem_layout, "TCBa", 4);

		GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_image_init: Allocated GMT Image %lx\n", (long)I);
		GMT_Report (API, GMT_MSG_DEBUG,
		            "gmtmex_image_init: Registered GMT Image array %lx via memory reference from MATLAB\n",
		            (long)I->data);
	}
	else {	/* Just allocate an empty container to hold an output image (signal this by passing 0s and NULLs) */
		if ((I = GMT_Create_Data (API, GMT_IS_IMAGE, GMT_IS_SURFACE, 0,
		                          NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_image_init: Failure to alloc GMT blank image container for holding output image\n");
#ifdef HAVE_GDAL
		GMT_Set_Default (API, "API_IMAGE_LAYOUT", "TCS");	/* State how we wish to receive images from GDAL */
#endif
	}
	return (I);
}

static void *gmtmex_dataset_init (void *API, unsigned int direction, unsigned int module_input, const mxArray *ptr) {
	/* Create containers to hold or receive data tables:
	 * direction == GMT_IN:  Create empty GMT_DATASET container, fill from Mex, and use as GMT input.
	 *	Input from MATLAB may be a MEX structure or a plain matrix
	 * direction == GMT_OUT: Create empty GMT_DATASET container, let GMT fill it out, and use for Mex output.
 	 * If direction is GMT_IN then we are given a MATLAB struct and can determine dimension.
	 * If output then we dont know size so we set dimensions to zero. */
	struct GMT_DATASET *D = NULL;

	if (direction == GMT_IN) {	/* Data given, dimensions are know, create container for GMT */
		uint64_t seg, col, start, k, n_headers, dim[4] = {1, 0, 0, 0};	/* We only return one table */
		size_t length = 0;
		char buffer[BUFSIZ] = {""};
		mxArray *mx_ptr = NULL;
		double *data = NULL;
		struct GMT_DATASEGMENT *S = NULL;

		if (!ptr) mexErrMsgTxt ("gmtmex_dataset_init: Input is empty where it can't be.\n");
		if (mxIsNumeric (ptr)) {	/* Got a MATLAB matrix as input - pass data pointers via MATRIX to save memory */
			struct GMT_MATRIX *M = NULL;
			unsigned int flag = (module_input) ? GMT_VIA_MODULE_INPUT : 0;
			mxClassID type = mxGetClassID (ptr);	/* Storage type for this matrix */
			dim[DIM_ROW] = mxGetM (ptr);		/* Number of rows */
			dim[DIM_COL] = mxGetN (ptr);		/* Number of columns */
			if ((M = GMT_Create_Data (API, GMT_IS_MATRIX|flag, GMT_IS_PLP, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
				mexErrMsgTxt ("gmtmex_dataset_init: Failure to alloc GMT source matrix\n");
			GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_dataset_init: Allocated GMT Matrix %lx\n", (long)M);
			switch (type) {	/* Assign ML type pointer to the corresponding GMT matrix union pointer */
				case mxDOUBLE_CLASS: M->type = GMT_DOUBLE; M->data.f8  =             mxGetData (ptr); break;
				case mxSINGLE_CLASS: M->type = GMT_FLOAT;  M->data.f4  =    (float *)mxGetData (ptr); break;
				case mxUINT64_CLASS: M->type = GMT_ULONG;  M->data.ui8 = (uint64_t *)mxGetData (ptr); break;
				case mxINT64_CLASS:  M->type = GMT_LONG;   M->data.si8 =  (int64_t *)mxGetData (ptr); break;
				case mxUINT32_CLASS: M->type = GMT_UINT;   M->data.ui4 = (uint32_t *)mxGetData (ptr); break;
				case mxINT32_CLASS:  M->type = GMT_INT;    M->data.si4 =  (int32_t *)mxGetData (ptr); break;
				case mxUINT16_CLASS: M->type = GMT_USHORT; M->data.ui2 = (uint16_t *)mxGetData (ptr); break;
				case mxINT16_CLASS:  M->type = GMT_SHORT;  M->data.si2 =  (int16_t *)mxGetData (ptr); break;
				case mxUINT8_CLASS:  M->type = GMT_UCHAR;  M->data.uc1 =  (uint8_t *)mxGetData (ptr); break;
				case mxINT8_CLASS:   M->type = GMT_CHAR;   M->data.sc1 =   (int8_t *)mxGetData (ptr); break;
				default:
					mexErrMsgTxt ("gmtmex_dataset_init: Unsupported MATLAB data type in GMT matrix input.");
					break;
			}
			/* Data from MATLAB and Octave(mex) is in col format and data from Octave(oct) is in row format */
#ifdef GMT_OCTOCT
			M->dim = M->n_columns;
#else
			M->dim = M->n_rows;
#endif
			M->alloc_mode = GMT_ALLOC_EXTERNALLY;	/* Since matrix was allocated by MATLAB/Octave we cannot free it in GMT */
			M->shape = MEX_COL_ORDER;		/* Either col or row order, depending on MATLAB/Octave setting in gmtmex.h */
			return (M);
		}
		/* We come here if we did not receive a matrix */
		if (!mxIsStruct (ptr)) mexErrMsgTxt ("gmtmex_dataset_init: Expected a data structure for input\n");
		dim[GMT_SEG] = mxGetM (ptr);	/* Number of segments */
		if (dim[GMT_SEG] == 0) mexErrMsgTxt ("gmtmex_dataset_init: Input has zero segments where it can't be.\n");
		if ((mx_ptr = mxGetField (ptr, 0, "data")) == NULL)	/* Get first segment's data matrix */
			mexErrMsgTxt("gmtmex_dataset_init: The 'data' array is NULL where it can't be\n");
		dim[GMT_COL] = mxGetN (mx_ptr);	/* Number of columns */
		if (dim[GMT_COL] == 0) mexErrMsgTxt ("gmtmex_dataset_init: Input has zero columns where it can't be.\n");
		if ((D = GMT_Create_Data (API, GMT_IS_DATASET, GMT_IS_PLP, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_dataset_init: Failure to alloc GMT destination dataset\n");
		GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_dataset_init: Allocated GMT dataset %lx\n", (long)D);

		for (seg = 0; seg < dim[GMT_SEG]; seg++) {	/* Each incoming structure is a new data segment */
			mx_ptr = mxGetField (ptr, (mwSize)seg, "header");	/* Get pointer to MEX segment header */
			buffer[0] = 0;	/* Reset our temporary text buffer */
			if ((length = mxGetN (mx_ptr)) != 0)	/* These is a non-empty segment header to keep */
				mxGetString (mx_ptr, buffer, (mwSize)(length+1));
			mx_ptr = mxGetField (ptr, (mwSize)seg, "data");	/* Data matrix for this segment */
			data = mxGetData (mx_ptr);
			dim[GMT_ROW] = mxGetM (mx_ptr);	/* Number of rows in matrix */
			/* Allocate a new data segment and hook up to table */
			S = GMT_Alloc_Segment (API, GMT_IS_DATASET, dim[GMT_ROW], dim[GMT_COL], buffer, D->table[0]->segment[seg]);
			for (col = start = 0; col < S->n_columns; col++, start += S->n_rows) /* Copy the data columns */
				memcpy (S->data[col], &data[start], S->n_rows * sizeof (double));
			D->table[0]->n_records += S->n_rows;	/* Must manually keep track of totals */
			if (seg == 0) {	/* First segment may have table information */
				mxArray *mx_ptr_t = mxGetField (ptr, (mwSize)seg, "comment");	/* Table headers */
				if (mx_ptr_t && (n_headers = mxGetM (mx_ptr_t)) != 0) {	/* Number of headers found */
					char *txt = NULL;
					for (k = 0; k < n_headers; k++) {
						mx_ptr = mxGetCell (mx_ptr_t, (mwSize)k);
						txt = mxArrayToString (mx_ptr);
						if (GMT_Set_Comment (API, GMT_IS_DATASET, GMT_COMMENT_IS_TEXT, txt, D))
							mexErrMsgTxt("gmtmex_dataset_init: Failed to set a dataset header\n");
					}
				}
			}
		}
		D->n_records = D->table[0]->n_records;
	}
	else {	/* Here we set up an empty container to receive data from GMT */
		if ((D = GMT_Create_Data (API, GMT_IS_DATASET, GMT_IS_PLP, 0, NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_dataset_init: Failure to alloc GMT source dataset\n");
		GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_dataset_init: Allocated GMT Dataset %lx\n", (long)D);
	}
	return (D);
}

static struct GMT_TEXTSET *gmtmex_textset_init (void *API, unsigned int direction, unsigned int module_input, unsigned int family, const mxArray *ptr) {
	/* Used to create containers to hold or receive mixed data/text tables:
	 * direction == GMT_IN:  Create empty GMT_TEXTSET container, fill, and use as GMT input.
	 *	input may be MEX segments, numerical matrix (will be converted to text), cell array of strings, or a single string.
	 * direction == GMT_OUT: Create empty GMT_TEXTSET container, let GMT fill it out, and use for Mex output.
 	 * If direction is GMT_IN then we are given a MATLAB struct and can determine dimension.
	 * If output then we dont know size so all we do is specify data type. */
	struct GMT_TEXTSET *T = NULL;
	
	if (direction == GMT_IN) {	/* Dimensions are known, extract them and set dim array for a GMT_MATRIX resource */
		uint64_t seg, col, row, n_cols = 0, k, n_headers, dim[3] = {1, 0, 0};	/* Only return one table */
		size_t length = 0;
		int add_text = 0;
		char buffer[BUFSIZ] = {""}, word[64] = {""}, *txt = NULL;
		mxArray *mx_ptr = NULL;
		double *data = NULL;
		struct GMT_TEXTSEGMENT *S = NULL;
		if (!ptr)
			mexErrMsgTxt("gmtmex_textset_init: Input pointer is NULL where it can't be.\n");
		else if (mxIsEmpty(ptr)) {
			mexPrintf("gmtmex_textset_init: Input text is empty, unknown consequence.\n");
			return NULL;
		}

		if (mxIsStruct (ptr)) {	/* Regular structure for data/textsets */
			mxArray *mx_ptr_d = NULL, *mx_ptr_t = NULL;
			dim[GMT_SEG] = mxGetM (ptr);	/* Number of segments */
			if (dim[GMT_SEG] == 0) mexErrMsgTxt ("gmtmex_textset_init: Input has zero segments where it can't be.\n");
			if ((T = GMT_Create_Data (API, GMT_IS_TEXTSET, GMT_IS_PLP, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
				mexErrMsgTxt ("gmtmex_textset_init: Failure to alloc GMT destination dataset\n");
			GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_textset_init: Allocated GMT textset %lx\n", (long)T);
			for (seg = 0; seg < dim[GMT_SEG]; seg++) {	/* Each incoming structure is a new segment */
				mx_ptr = mxGetField (ptr, (mwSize)seg, "header");	/* Segment header */
				buffer[0] = 0;	/* Reset our temporary buffer */
				if ((length = mxGetN (mx_ptr)) != 0)
					mxGetString (mx_ptr, buffer, (mwSize)(length+1));
				mx_ptr_d = mxGetField (ptr, (mwSize)seg, "data");	/* Data table for this segment */
				mx_ptr_t = mxGetField (ptr, (mwSize)seg, "text");	/* Text table for this segment */
				if (mxIsEmpty(mx_ptr_t))	/* The cell array is empty so no text will be converted */
					add_text = 0;
				else {	/* Get pointer to cell array and get row count */
					dim[GMT_ROW] = mxGetM (mx_ptr_t);	/* Number of rows found */
					add_text = 1;
				}
				if (dim[GMT_ROW] == 0)	/* No text array present, rely on data array instead. */
					dim[GMT_ROW] = mxGetM (mx_ptr_d);	/* Number of rows */
				n_cols = mxGetN (mx_ptr_d);	/* Number of data cols, if any */
				/* Allocate new text segment and hook it up to the table */
				S = GMT_Alloc_Segment (API, GMT_IS_TEXTSET, dim[GMT_ROW], 0, buffer, T->table[0]->segment[seg]);
				data = mxGetData (mx_ptr_d);
				/* Combine any data and cell arrays into text records */
				for (row = 0; row < S->n_rows; row++) {
					/* First deal with the [optional] data matrix for leading columns */
					if (n_cols) sprintf (buffer, "%.16g", data[row]);
					for (col = 1; col < n_cols; col++) {
						sprintf (word, "\t%.16g", data[row+col*S->n_rows]);
						strcat (buffer, word);
					}
					if (add_text) {	/* Then append the optional text strings */
						mx_ptr = mxGetCell (mx_ptr_t, (mwSize)row);	/* Pointer to this cell */
						if ((txt = mxArrayToString (mx_ptr)) != NULL) {		/* Yes, there was a string there */
							strcat (buffer, "\t");	/* Append the string after a tab */
							strcat (buffer, txt);
							mxFree (txt);
						}
					}
					S->data[row] = GMT_Duplicate_String (API, buffer);
				}
				if (seg == 0) {	/* First segment may have dataset information */
					mx_ptr_t = mxGetField (ptr, (mwSize)seg, "comment");	/* Table headers */
					if (mx_ptr_t && (n_headers = mxGetM (mx_ptr_t)) != 0) {	/* Number of headers found */
						for (k = 0; k < n_headers; k++) {	/* Add each header as a text comment to the textset */
							mx_ptr = mxGetCell (mx_ptr_t, (mwSize)row);
							txt = mxArrayToString (mx_ptr);
							if (GMT_Set_Comment (API, GMT_IS_TEXTSET, GMT_COMMENT_IS_TEXT, txt, T))
								mexErrMsgTxt("gmtmex_textset_init: Failed to set a textset header\n");
						}
					}
				}
			}
		}
		else {	/* Get here when given either a numerical matrix, a cell array or a single string */
			unsigned int got_text = 0;
			dim[GMT_ROW] = mxGetM(ptr); /* Number of records */
			dim[GMT_SEG] = 1;           /* Only one segment is returned here */
			if (mxIsNumeric (ptr)) {    /* Got data matrix instead, must convert to text first (below) */
				n_cols = mxGetN (ptr);  /* Number of columns */
				data = mxGetData (ptr); /* Get pointer to the data matrix */
			}
			else if (mxIsChar(ptr) && dim[GMT_ROW] == 1)	/* Special case: Got a single text record instead of cell array */
				got_text = 1;
			else if (!mxIsCell (ptr))	/* Merda */
				mexErrMsgTxt ("gmtmex_textset_init: Expected either a Cell array, Matrix, or text string for input\n");
			if (n_cols == 0 && dim[GMT_ROW] == 1 && !got_text) {	/* Check if we got a transpose arrangement or just one record */
				row = mxGetN (ptr);                 /* Also possibly number of records */
				if (row > 1) dim[GMT_ROW] = row;    /* User gave row-vector of cells */
			}
			if ((T = GMT_Create_Data (API, family, GMT_IS_NONE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
				mexErrMsgTxt ("gmtmex_textset_init: Failure to alloc GMT source TEXTSET for input\n");
			S = T->table[0]->segment[0];	/* Only one segment (already allocated and hooked up) will be returned by MATLAB */
			S->n_rows = dim[GMT_ROW];
			T->alloc_mode = GMT_ALLOC_EXTERNALLY;
			if (got_text)	/* Just got that single record */
				S->data[0] = mxArrayToString (ptr);
			else {	/* Must get strings out of the cell array or reformat data matrix */
				for (row = 0; row < S->n_rows; row++) {
					if (n_cols) {	/* Create text string from matrix */
						sprintf (buffer, "%.16g", data[row]);
						for (col = 1; col < n_cols; col++) {
							sprintf (word, "\t%.16g", data[row+col*dim[GMT_ROW]]);
							strcat (buffer, word);
						}
						txt = mxstrdup (buffer);
					}
					else {	/* Got cell array */
						mx_ptr = mxGetCell (ptr, (mwSize)row);
						txt = mxArrayToString (mx_ptr);
					}
					S->data[row] = txt;
				}
			}
			T->n_records = T->table[0]->n_records = S->n_rows;
		}
		GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_textset_init: Allocated GMT TEXTSET %lx\n", (long)T);
	}
	else {	/* Here we set up an empty container to receive a textset from GMT */
		if ((T = GMT_Create_Data (API, GMT_IS_TEXTSET, GMT_IS_NONE, 0, NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_textset_init: Failure to alloc GMT source textset\n");
		GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_textset_init: Allocated GMT Textset %lx\n", (long)T);
	}
	return (T);
}

static struct GMT_PALETTE *gmtmex_cpt_init (void *API, unsigned int direction, unsigned int module_input, const mxArray *ptr) {
	/* Used to create an empty CPT container to hold a GMT Color Palette.
 	 * If direction is GMT_IN then we are given a MATLAB CPT struct and can determine its size, etc.
	 * If direction is GMT_OUT then we allocate an empty GMT CPT as a destination. */
	struct GMT_PALETTE *P = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the input pointer */
		unsigned int k, j, one = 1, n_headers, *depth = NULL;
		uint64_t dim[2] = {0, 0};
		unsigned int flag = (module_input) ? GMT_VIA_MODULE_INPUT : 0;
		char model[8] = {""};
		mxArray *mx_ptr[N_MEX_FIELDNAMES_CPT];
		double *colormap = NULL, *range = NULL, *minmax = NULL, *alpha = NULL, *bfn = NULL, *hinge = NULL, *cpt = NULL;

		if (mxIsEmpty (ptr))
			mexErrMsgTxt ("gmtmex_cpt_init: The input that was supposed to contain the CPT, is empty\n");
		if (!mxIsStruct (ptr))
			mexErrMsgTxt ("gmtmex_cpt_init: Expected a CPT structure for input\n");
		for (k = 0; k < N_MEX_FIELDNAMES_CPT; k++) {
			if ((mx_ptr[k] = mxGetField (ptr, 0, GMTMEX_fieldname_cpt[k])) == NULL)
				gmtmex_quit_if_missing ("gmtmex_cpt_init", GMTMEX_fieldname_cpt[k]);
		}
		
		dim[0] = mxGetM (mx_ptr[0]);	/* Number of rows in colormap */
		if (dim[0] < 1)
			mexErrMsgTxt ("gmtmex_cpt_init: Colormap array has no CPT values\n");
		colormap = mxGetData (mx_ptr[0]);
		alpha    = mxGetData (mx_ptr[1]);
		range    = mxGetData (mx_ptr[2]);
		minmax   = mxGetData (mx_ptr[3]);
		bfn      = mxGetData (mx_ptr[4]);
		depth    = mxGetData (mx_ptr[5]);
		hinge    = mxGetData (mx_ptr[6]);
		cpt      = mxGetData (mx_ptr[7]);

		dim[1] = mxGetM (mx_ptr[2]);	/* Length of range array */
		if (dim[0] > dim[1]) {	/* This only happens when we have a continuous color table */
			dim[1] = dim[0];    /* Actual length of colormap array */
			dim[0]--;           /* Number of CPT slices */
		}
		else	/* Discrete, so the one offset needs to be zero */
			one = 0;

		if ((P = GMT_Create_Data (API, GMT_IS_PALETTE|flag, GMT_IS_NONE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_cpt_init: Failure to alloc GMT source CPT for input\n");

		if ((n_headers = (unsigned int)mxGetM (mx_ptr[9])) != 0) {	/* Number of headers found */
			char *txt = NULL;
			mxArray *ptr = NULL;
			for (k = 0; k < n_headers; k++) {
				ptr = mxGetCell (mx_ptr[9], (mwSize)k);
				txt = mxArrayToString (ptr);
				if (GMT_Set_Comment (API, GMT_IS_PALETTE, GMT_COMMENT_IS_TEXT, txt, P))
					mexErrMsgTxt("gmtmex_cpt_init: Failed to set a CPT header\n");
			}
		}
		for (j = 0; j < 3; j++) {	/* Do the bfn first */
			for (k = 0; k < 3; k++)
				P->bfn[j].rgb[k] = bfn[j+k*3];
		}
		for (j = 0; j < P->n_colors; j++) {	/* OK to access j+1'th elemenent since length of colormap is P->n_colors+1 */
			for (k = 0; k < 3; k++) {
				P->data[j].rgb_low[k]  = cpt[j+k*dim[0]];
				P->data[j].rgb_high[k] = cpt[j+(k+3)*dim[0]];
			}
			P->data[j].rgb_low[3]  = alpha[j];
			P->data[j].rgb_high[3] = alpha[j+one];
			P->data[j].z_low  = range[j];
			P->data[j].z_high = range[j+P->n_colors];
			P->data[j].annot = 3;	/* Enforce annotations for now */
		}
		P->is_continuous = one;
		P->is_bw = P->is_gray = 0;
		if (depth[0] == 1)
			P->is_bw = 1;
		else if (depth[0] == 8)
			P->is_gray = 1;
		GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_cpt_init: Allocated GMT CPT %lx\n", (long)P);
		if (mxIsNaN (hinge[0])) {
			P->has_hinge = 1;
			P->mode &= GMT_CPT_HINGED;
		}
		mxGetString (mx_ptr[8], model, (mwSize)mxGetN(mx_ptr[8])+1);
		if (!strncmp (model, "hsv", 3U))
			P->model = GMT_HSV;
		else if (!strncmp (model, "cmyk", 4U))
			P->model = GMT_CMYK;
		else
			P->model = GMT_RGB;
	}
	else {	/* Just allocate an empty container to hold an output grid (signal this by passing 0s and NULLs) */
		if ((P = GMT_Create_Data (API, GMT_IS_PALETTE, GMT_IS_NONE, 0, NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_cpt_init: Failure to alloc GMT blank CPT container for holding output CPT\n");
	}
	return (P);
}

#if GMT_MINOR_VERSION > 2
static struct GMT_POSTSCRIPT *gmtmex_ps_init (void *API, unsigned int direction, unsigned int module_input, const mxArray *ptr) {
	/* Used to Create an empty POSTSCRIPT container to hold a GMT POSTSCRIPT object.
 	 * If direction is GMT_IN then we are given a MATLAB structure with known sizes.
	 * If direction is GMT_OUT then we allocate an empty GMT POSTSCRIPT as a destination. */
	struct GMT_POSTSCRIPT *P = NULL;
	if (direction == GMT_IN) {	/* Dimensions are known from the MATLAB input pointer */
		uint64_t dim[1] = {0}, *length = NULL;
		unsigned int k, n_headers, *mode = NULL, flag = (module_input) ? GMT_VIA_MODULE_INPUT : 0;
		mxArray *mx_ptr[N_MEX_FIELDNAMES_PS];
		char *PS = NULL;

		if (mxIsEmpty (ptr))
			mexErrMsgTxt ("gmtmex_ps_init: The input that was supposed to contain the PostScript structure is empty\n");
		if (!mxIsStruct (ptr))
			mexErrMsgTxt ("gmtmex_ps_init: Expected a MATLAB PostScript structure for input\n");
		for (k = 0; k < N_MEX_FIELDNAMES_PS; k++) {
			if ((mx_ptr[k] = mxGetField (ptr, 0, GMTMEX_fieldname_ps[k])) == NULL)
				gmtmex_quit_if_missing ("gmtmex_ps_init", GMTMEX_fieldname_ps[k]);
		}
			
		length = mxGetData (mx_ptr[1]);
		if (length[0] == 0)
			mexErrMsgTxt ("gmtmex_ps_init: Dimension of PostScript given as zero\n");
		PS = malloc (mxGetN(mx_ptr[0])+1);
		mxGetString (mx_ptr[0], PS, (mwSize)mxGetN(mx_ptr[0]));
		mode = mxGetData (mx_ptr[2]);
		/* Passing dim[0] = 0 since we dont want any allocation of a PS string */
		if ((P = GMT_Create_Data (API, GMT_IS_POSTSCRIPT|flag, GMT_IS_NONE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_ps_init: Failure to alloc GMT POSTSCRIPT source for input\n");
		P->data = PS;	/* PostScript string instead is coming from MATLAB */
		P->alloc_mode = GMT_ALLOC_EXTERNALLY;	/* Hence we are not allowed to free it */
		P->n_bytes = length[0];	/* Length of the actual PS string */
		P->n_alloc = 0;		/* But nothing was actually allocated here - just passing pointer from MATLAB */
		P->mode = mode[0];	/* Inherit the mode */
		if ((n_headers = (unsigned int)mxGetM (mx_ptr[3])) != 0) {	/* Number of headers found */
			char *txt = NULL;
			mxArray *ptr = NULL;
			for (k = 0; k < n_headers; k++) {
				ptr = mxGetCell (mx_ptr[3], (mwSize)k);
				txt = mxArrayToString (ptr);
				if (GMT_Set_Comment (API, GMT_IS_POSTSCRIPT, GMT_COMMENT_IS_TEXT, txt, P))
					mexErrMsgTxt("gmtmex_ps_init: Failed to set a PostScript header\n");
			}
		}
		GMT_Report (API, GMT_MSG_DEBUG, "gmtmex_ps_init: Allocated GMT POSTSCRIPT %lx\n", (long)P);
	}
	else {	/* Just allocate an empty container to hold an output PS object (signal this by passing 0s and NULLs) */
		if ((P = GMT_Create_Data (API, GMT_IS_POSTSCRIPT, GMT_IS_NONE, 0, NULL, NULL, NULL, 0, 0, NULL)) == NULL)
			mexErrMsgTxt ("gmtmex_ps_init: Failure to alloc GMT POSTSCRIPT container for holding output PostScript\n");
	}
	return (P);
}
#endif

void GMTMEX_objecttype (char *type, const mxArray *ptr) {
	/* Determine what we are returning so gmt write can pass the correct -T? flag */
	mxArray *mx_ptr = NULL;
	if (mxIsEmpty (ptr))
		mexErrMsgTxt ("GMTMEX_objecttype: Pointer is empty\n");
	if (mxIsStruct (ptr)) {	/* So either dataset, grid, image, cpt, or PS */
		mx_ptr = mxGetField (ptr, 0, "data");
		if (mx_ptr) { type[3] = 'd'; return;}
		mx_ptr = mxGetField (ptr, 0, "postscript");
		if (mx_ptr) { type[3] = 'p'; return;}
		mx_ptr = mxGetField (ptr, 0, "hinge");
		if (mx_ptr) {type[3] = 'c'; return;}
		mx_ptr = mxGetField (ptr, 0, "image");
		if (mx_ptr) {type[3] = 'i'; return;}
		mx_ptr = mxGetField (ptr, 0, "z");
		if (mx_ptr) {type[3] = 'g'; return;}
		mexErrMsgTxt ("GMTMEX_objecttype: Could not recognize the structure\n");
	}
	else if (mxIsCell (ptr))
		type[3] = 't';
	else
		type[3] = 'd';
}

void *GMTMEX_Register_IO (void *API, struct GMT_RESOURCE *X, const mxArray *ptr) {
	/* Create the grid or matrix container, register it, and return the ID */
	void *obj = NULL;	/* Pointer to the container we created */
	char *name[2] = {"Matrix", "CellArray"};
	unsigned int module_input = (X->option->option == GMT_OPT_INFILE);
	X->object_ID = GMT_NOTSET;

	switch (X->family) {
		case GMT_IS_GRID:	/* Get a grid from Matlab or a dummy one to hold GMT output */
			obj = gmtmex_grid_init (API, X->direction, module_input, ptr);
			X->object_ID = GMT_Get_ID (API, GMT_IS_GRID, X->direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got Grid with Object ID %d\n", X->object_ID);
			break;
		case GMT_IS_IMAGE:	/* Get an image from Matlab or a dummy one to hold GMT output */
			obj = gmtmex_image_init (API, X->direction, module_input, ptr);
			X->object_ID = GMT_Get_ID (API, GMT_IS_IMAGE, X->direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got Image with Object ID %d\n", X->object_ID);
			break;
		case GMT_IS_DATASET:	/* Get a dataset from Matlab or a dummy one to hold GMT output */
			/* Ostensibly a DATASET, but it might be a TEXTSET passed via a cell array, so we must check */
			if (X->direction == GMT_IN && (mxIsCell (ptr) || mxIsChar (ptr))) {	/* Got text input */
				obj = gmtmex_textset_init (API, X->direction, module_input, GMT_IS_TEXTSET, ptr);
				X->family = GMT_IS_TEXTSET;
			}
			else	/* Got something for which a dataset container is appropriate */
				obj = gmtmex_dataset_init (API, X->direction, module_input, ptr);
			X->object_ID = GMT_Get_ID (API, X->family, X->direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got %s with Object ID %d\n", name[X->family], X->object_ID);
			break;
		case GMT_IS_TEXTSET:	/* Get a textset from Matlab or a dummy one to hold GMT output */
			obj = gmtmex_textset_init (API, X->direction, module_input, GMT_IS_TEXTSET, ptr);
			X->object_ID = GMT_Get_ID (API, GMT_IS_TEXTSET, X->direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got TEXTSET with Object ID %d\n", X->object_ID);
			break;
		case GMT_IS_PALETTE:	/* Get a palette from Matlab or a dummy one to hold GMT output */
			obj = gmtmex_cpt_init (API, X->direction, module_input, ptr);
			X->object_ID = GMT_Get_ID (API, GMT_IS_PALETTE, X->direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got CPT with Object ID %d\n", X->object_ID);
			break;
#if GMT_MINOR_VERSION > 2
		case GMT_IS_POSTSCRIPT:	/* Get a PostScript struct from Matlab or a dummy one to hold GMT output */
			obj = gmtmex_ps_init (API, X->direction, module_input, ptr);
			X->object_ID = GMT_Get_ID (API, GMT_IS_POSTSCRIPT, X->direction, obj);
			GMT_Report (API, GMT_MSG_DEBUG, "GMTMEX_Register_IO: Got POSTSCRIPT with Object ID %d\n", X->object_ID);
			break;
#endif
		default:
			GMT_Report (API, GMT_MSG_NORMAL, "GMTMEX_Register_IO: Bad data type (%d)\n", X->family);
			break;
	}
	return (obj);
}
