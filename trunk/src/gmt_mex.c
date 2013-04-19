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
#include <string.h>
#include <ctype.h>
#include <math.h>

int GMTMEX_print_func (FILE *fp, const char *message)
{
	/* Repalcement for GMT's gmt_print_func.  It is being used indirectly via
	 * API->print_func.  Purpose of this is to allow Matlab (which cannot use
	 * printf) to reset API->print_func this functions via GMT_Create_Session. */

	mexErrMsgTxt (message);
	return 0;
}

void GMTMEX_grdheader2info (mxArray *plhs[], struct GMT_GRID *G, int item)
{	/* Return the grid's header info via the info array */
	double *hdr = NULL;
	plhs[item] = mxCreateDoubleMatrix (1, 9, mxREAL);
	hdr = mxGetPr (plhs[item]);
	memcpy (hdr, G->header->wesn, 4 * sizeof(double));
	hdr[4] = G->header->z_min;	hdr[5] = G->header->z_max;
	hdr[6] = (double)G->header->registration;
	hdr[7] = G->header->inc[GMT_X];	hdr[8] = G->header->inc[GMT_Y];
}

double *GMTMEX_info2grdheader (void *API, const mxArray *prhs[], int nrhs, struct GMT_GRID_HEADER *header)
{	/* Return the grid's header info via the info array (nrhs == 3) or x/y arrays (nrhs == 4|5) */
	double *z = NULL;
	if (nrhs == 3) {	/* Gave Z, info */
		double *hdr = NULL;
		z = mxGetData (prhs[0]);
		header->nx = mxGetN (prhs[0]);
		header->ny = mxGetM (prhs[0]);
		hdr = mxGetData (prhs[1]);
		memcpy (header->wesn, hdr, 4 * sizeof(double));
		header->z_min = hdr[4];
		header->z_max = hdr[5];
		header->registration = lrint (hdr[6]);
		header->inc[GMT_X] = hdr[7];
		header->inc[GMT_Y] = hdr[8];
	}
	else {	/* Gave x, y, Z [reg] */
		double *r = NULL, *x = NULL, *y = NULL;
		int col, row, error = 0;
		x = mxGetData (prhs[0]);
		y = mxGetData (prhs[1]);
		z = mxGetData (prhs[2]);
		header->nx = mxGetN (prhs[2]);
		header->ny = mxGetM (prhs[2]);
		header->inc[GMT_X] = x[1] - x[0];
		header->inc[GMT_Y] = y[1] - y[0];
		for (col = 2; !error && col < header->nx; col++) 
			if ((x[col] - x[col-1]) != header->inc[GMT_X]) error = 1;
		for (row = 2; !error && row < header->ny; row++) 
			if ((y[row] - y[row-1]) != header->inc[GMT_Y]) error = 1;
		if (error) {
			mexErrMsgTxt ("grdwrite: x and/or y not equidistant");
		}
		if (nrhs == 5) {
			r = mxGetData (prhs[3]);
			header->registration = lrint (r[0]);
		}
		else
			header->registration = GMT_GRID_NODE_REG;
		header->wesn[GMT_XLO] = (header->registration == GMT_GRID_PIXEL_REG) ? x[0] - 0.5 * header->inc[GMT_X] : x[0];
		header->wesn[GMT_XHI] = (header->registration == GMT_GRID_PIXEL_REG) ? x[header->nx-1] + 0.5 * header->inc[GMT_X] : x[header->nx-1];
		header->wesn[GMT_YLO] = (header->registration == GMT_GRID_PIXEL_REG) ? y[0] - 0.5 * header->inc[GMT_Y] : y[0];
		header->wesn[GMT_YHI] = (header->registration == GMT_GRID_PIXEL_REG) ? y[header->ny-1] + 0.5 * header->inc[GMT_Y] : y[header->ny-1];
	}
	return (z);
}

void GMTMEX_grdxy (void *API, mxArray *plhs[], struct GMT_GRID *G, int px, int py)
{	/* Return x,y arrays also */
	unsigned int row;
	double *xg = NULL, *yg = NULL, *tmp = NULL;
	plhs[px] = mxCreateDoubleMatrix (1, G->header->nx, mxREAL);
	plhs[py] = mxCreateDoubleMatrix (1, G->header->ny, mxREAL);
	xg = mxGetPr (plhs[px]);	yg = mxGetPr (plhs[py]);
	/* Fill in the x and y arrays; note Matlab y array is flipped relative to the GMT y array */
	tmp = GMT_Get_Coord (API, GMT_IS_GRID, GMT_X, G);
	memcpy (xg, tmp, G->header->nx * sizeof (double));
	GMT_Destroy_Data (API, GMT_CLOBBER, &tmp);
	tmp = GMT_Get_Coord (API, GMT_IS_GRID, GMT_Y, G);
	for (row = 0; row < G->header->ny; row++) yg[G->header->ny-row-1] = tmp[row];
	GMT_Destroy_Data (API, GMT_CLOBBER, &tmp);
}

char *GMTMEX_src_vector_init (void *API, const mxArray *prhs[], unsigned int n_cols, int n_start, struct GMT_VECTOR **V)
{	/* Used by programs that expect either an input file name or data vectors x, y[, other cols] */
	char *i_string = NULL;
	if (mxIsChar(prhs[0]))		/* Gave a file name */
		i_string = mxArrayToString (prhs[0]);	/* Load the file name into a char string */
 	else {				/* Input via two or more column vectors */
		int col, in_ID;
		uint64_t dim[1] = {n_cols};
		i_string = mxMalloc (BUFSIZ);
		if ((*V = GMT_Create_Data (API, GMT_IS_VECTOR, GMT_IS_NONE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL) mexErrMsgTxt ("Failure to alloc GMT source vectors\n");
		for (col = n_start; col < n_cols+n_start; col++) {	/* Hook up one vector per column and determine data type */
			if (mxIsDouble(prhs[col])) {
				(*V)->type[col] = GMT_DOUBLE;
				(*V)->data[col].f8 = mxGetData (prhs[col]);
			}
			else if (mxIsSingle(prhs[col])) {
				(*V)->type[col] = GMT_FLOAT;
				(*V)->data[col].f4 = (float *)mxGetData (prhs[col]);
			}
			else if (mxIsInt32(prhs[col])) {
				(*V)->type[col] = GMT_INT;
				(*V)->data[col].si4 = (int32_t *)mxGetData (prhs[col]);
			}
			else if (mxIsInt16(prhs[col])) {
				(*V)->type[col] = GMT_SHORT;
				(*V)->data[col].si2 = (int16_t *)mxGetData (prhs[col]);
			}
			else if (mxIsInt8(prhs[col])) {
				(*V)->type[col] = GMT_CHAR;
				(*V)->data[col].sc1 = (int8_t *)mxGetData (prhs[col]);
			}
			else
				mexErrMsgTxt ("Unsupported data type in GMT input.");
		}

		(*V)->n_rows = MAX (mxGetM (prhs[0]), mxGetN (prhs[0]));	/* So it works for both column or row vectors */
		(*V)->n_columns = n_cols;
		if ((in_ID = GMT_Register_IO (API, GMT_IS_DATASET, GMT_IS_READONLY + GMT_VIA_VECTOR, GMT_IS_NONE, GMT_IN, NULL, V)) == GMTAPI_NOTSET) {
			mexErrMsgTxt ("Failure to register GMT source vectors\n");
		}
		if (GMT_Encode_ID (API, i_string, in_ID) != GMT_NOERROR) {		/* Make filename with embedded object ID */
			mexErrMsgTxt ("GMTMEX_parser: Failure to encode string\n");
		}
	}
	return (i_string);
}

char *GMTMEX_src_grid_init (void *API, const mxArray *prhs[], int nrhs)
{	/* Used by programs that expect either an input file name or matrix with info or data vectors x, y */
	char *i_string = NULL;
	struct GMT_GRID *G = NULL;
	if (nrhs == 2)		/* Gave a file name */
		i_string = mxArrayToString (prhs[0]);	/* Load the file name into a char string */
 	else {			/* Input via matrix and either info array or x,y arrays */
		unsigned int row, col;
		int in_ID;
		uint64_t gmt_ij;
		double *z = NULL;
		struct GMT_GRID_HEADER h_tmp;
		i_string = mxMalloc (BUFSIZ);

		/*  Get the Z array and fill in the header info */
		z = GMTMEX_info2grdheader (API, prhs, nrhs, &h_tmp);

		if ((G = GMT_Create_Data (API, GMT_IS_GRID, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, h_tmp.wesn, h_tmp.inc, \
			h_tmp.registration, GMTAPI_NOTSET, NULL)) == NULL) mexErrMsgTxt ("Failure to alloc GMT grid\n");
		
		/* Transpose from Matlab orientation to grd orientation */
		for (gmt_ij = row = 0; row < h_tmp.ny; row++) for (col = 0; col < h_tmp.nx; col++, gmt_ij++)
			G->data[gmt_ij] = (float)z[MEX_IJ(G,row,col)];
		if ((in_ID = GMT_Register_IO (API, GMT_IS_GRID, GMT_IS_REFERENCE, GMT_IS_SURFACE, GMT_IN, NULL, G)) == GMTAPI_NOTSET) {
			mexErrMsgTxt ("Failure to register GMT source grid\n");
		}
		if (GMT_Encode_ID (API, i_string, in_ID) != GMT_NOERROR) {	/* Make filename with embedded object ID */
			mexErrMsgTxt ("GMTMEX_parser: Failure to encode string\n");
		}
	}
	return (i_string);
}

char *GMTMEX_dest_grid_init (void *API, int *out_ID, int nlhs, char *options)
{	/* Associate output grid with Matlab grid */
	char *o_string = NULL;
	if ((*out_ID = GMT_Register_IO (API, GMT_IS_GRID, GMT_IS_REFERENCE, GMT_IS_SURFACE, GMT_OUT, NULL, NULL)) == GMTAPI_NOTSET) {
		mexErrMsgTxt ("Failure to register GMT destination grid\n");
	}
	if (nlhs == 0) {
		if (strstr (options, "-G")) 	/* User gave -G<file> among the options */
			return (NULL);		/* No output will be send to Matlab */
		else
			mexErrMsgTxt ("Error: neither -G option nor left hand side output args.");
	}
	o_string = mxMalloc(GMTAPI_STRLEN);
	if (GMT_Encode_ID (API, o_string, *out_ID) != GMT_NOERROR) {	/* Make filename with embedded object ID */
		mexErrMsgTxt ("GMTMEX_parser: Failure to encode string\n");
	}
	return (o_string);
}

char *GMTMEX_dest_vector_init (void *API, unsigned int n_cols, struct GMT_VECTOR **V, int nlhs, char *options)
{	/* Associate output data with Matlab/Octave vectors */
	char *o_string = NULL;
	unsigned int col;
	int out_ID;
	uint64_t dim[1] = {n_cols};

	o_string = mxMalloc(GMTAPI_STRLEN);
	if (nlhs == 0) {
		if (strstr (options, ">")) 	/* User gave > file among the options */
			return (NULL);		/* No output will be send to Matlab */
		else
			mexErrMsgTxt ("Error: neither output file name with the '>' "
					"redirection operator nor left hand side output args.");
	}
	if ((*V = GMT_Create_Data (API, GMT_IS_VECTOR, GMT_IS_NONE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL) mexErrMsgTxt ("Failure to alloc GMT source vectors\n");
	for (col = 0; col < n_cols; col++) (*V)->type[col] = GMT_DOUBLE;
	(*V)->alloc_mode = GMT_REFERENCE;
	if ((out_ID = GMT_Register_IO (API, GMT_IS_DATASET, GMT_IS_REFERENCE + GMT_VIA_VECTOR, GMT_IS_NONE, GMT_OUT, NULL, V)) == GMTAPI_NOTSET) {
		mexErrMsgTxt ("Failure to register GMT destination vectors\n");
	}
		
	if (GMT_Encode_ID (API, o_string, out_ID) != GMT_NOERROR) {	/* Make filename with embedded object ID */
		mexErrMsgTxt ("GMTMEX_parser: Failure to encode string\n");
	}
	return (o_string);
}

void GMTMEX_prep_mexgrd (void *API, mxArray *plhs[], int nlhs, struct GMT_GRID *G)
{	/* Turn a GMT grid into a 2-D Z matrix (and possibly x,y arrays) for passing back to Matlab/Octave */
	int px = -1, py = -1, pz = -1, pi = -1;
	unsigned int row, col;
	uint64_t gmt_ij;
	float *z = NULL;
	
	pz = (nlhs >= 3) ? 2 : 0;
	if (nlhs == 2 || nlhs == 4) pi = nlhs - 1;
	if (nlhs > 2) {px = 0; py = 1;}

	/* A. Create 2-D matrices for the return matrix */
	plhs[pz] = mxCreateNumericMatrix (G->header->ny, G->header->nx, mxSINGLE_CLASS, mxREAL);
	z = mxGetData (plhs[pz]);

	/* B. Load the real grd array into a double matlab array by
              transposing from padded GMT grd format to unpadded matlab format */

	for (gmt_ij = row = 0; row < G->header->ny; row++) for (col = 0; col < G->header->nx; col++, gmt_ij++)
		z[MEX_IJ(G,row,col)] = G->data[gmt_ij];
    
	/* C. Create header and x,y arrays, if requested  */

	if (pi >= 0) GMTMEX_grdheader2info (plhs, G, pi);	/* Also return info array */
	if (px >= 0) GMTMEX_grdxy (API, plhs, G, px, py);	/* Return x,y arrays also */
}

void GMTMEX_prep_mextbl (void *API, mxArray *plhs[], int nlhs, struct GMT_VECTOR *V)
{	/* We must duplicate output vectors since Matlab won't allow mixing of allocated arrays from outside Matlab */
	int p;
	double *z = NULL;
	for (p = 0; p < nlhs; p++) {
		plhs[p] = mxCreateNumericMatrix (V->n_rows, 1, mxDOUBLE_CLASS, mxREAL);
		z = mxGetData (plhs[p]);
		memcpy (z, V->data[p].f8, V->n_rows * sizeof(double));
	}
}

#if 0
char *GMTMEX_options_init (void *API, const mxArray *prhs[], int nrhs)
{	/* Secure string with the user's program options and parse -V, if included */
	int k = 1, j;
	char *options = NULL, *s = NULL;

	if (!mxIsChar(prhs[nrhs-1])) return (NULL);	/* No options in this case */
	options = mxArrayToString (prhs[nrhs-1]);

	if ((s = strstr (options, "-V"))) {	/* User gave -V[level] among the options */
		int level;
		s += 2;	/* Skip to char after -V */
		level = (isdigit ((int)s[0]) && s[0] >= '0' && s[0] <= '4') ? (int)(s[0] - '0') : GMT_MSG_VERBOSE;
		API->GMT->current.setting.verbose = level;
	}
	while (options[k] && !(options[k] == '>' && options[k-1] == ' ')) k++;		/* Test if ... */
	if (k < strlen(options) - 1) {			/* User used the "... > file ..." construct */
		int	len;
		len = strlen(options);
		options = mxRealloc(options, len+1);	/* Increase by one because we need to insert a '-' */
		for (j = len; j > k; j--)
			options[j] = options[j-1];	/* Make room to the to-be-inseted '-' */
		k++;
		options[k-1] = '-';
		while (options[k+1] == ' ') {	/* Remove all spaces between the '>' and filename */
			for (j = k+1; j < strlen(options)-1; j++)
				options[j] = options[j+1];
			options[j] = '\0';
		}
	}
	/* Now test for this mistake "... -> file ...". Note that now no need to expand the options string */
	else if ((s = strstr (options, "-> "))) {
		while (s[2] == ' ') {		/* Remove all spaces between the '>' and filename */
			for (j = 2; j < strlen(s)-1; j++)
				s[j] = s[j+1];
			s[j] = '\0';
		}
	}
	return (options);
}
#endif
#if 0
char *GMTMEX_build_cmd (void *API, char *src, char *options, char *dest, int mode)
{	/* Create the command based on options, src, and dist, which depends slightly on output type */
	char *cmd = mxMalloc (BUFSIZ);
	if (mode == GMT_IS_GRID) {
		if (dest)
			sprintf (cmd, "%s %s -G%s", src, options, dest);
		else
			sprintf (cmd, "%s %s", src, options);
	}
	else if (mode == GMT_IS_DATASET)
		if (dest)
			sprintf (cmd, "%s %s > %s", src, options, dest);
		else
			sprintf (cmd, "%s %s", src, options);
	else {	/* PostScript */
		(src) ? sprintf (cmd, "%s %s", src, options) : sprintf (cmd, "%s", options);
	}
	return (cmd);
}
#endif

void GMTMEX_free (char *input, char *output, char *options, char *cmd) {
	/* Free temporary local variables */
	if (input) mxFree (input);
	if (output) mxFree (output);	
	if (options) mxFree (options);	
	mxFree (cmd);
}
