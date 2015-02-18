/*
 *	$Id$
 *
 *	Copyright (c) 1991-2015 by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis and F. Wobbe
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
/* GMT convenience functions used by MATLAB/OCTAVE mex/oct functions
 */

#ifndef GMTMEX_H
#define GMTMEX_H

#include "gmt.h"
#include <string.h>
#include <ctype.h>

#ifdef NO_MEX	/* THis would just be for testing the parser */
#define mxArray void
char revised_cmd[BUFSIZ];	/* Global variable used to show revised command when testing only */
#else
#ifdef GMT_OCTOCT
#include <oct.h>
#else
#include <mex.h>
#define mxIsScalar_(mx) \
	( (2 == mxGetNumberOfDimensions(mx)) \
		&&  (1 == mxGetM(mx))&&  (1 == mxGetN(mx)) )
#endif	/* Matlab and Octave(mex) */
#endif	/* NO_MEX */

/* Matlab and Octave (in -mex mode) are identical, oct files are different and not yet tested */

#ifdef GMT_OCTOCT	/* Octave oct files only */
#define MEX_PROG "Octave(oct)"
#define MEX_COL_ORDER GMT_IS_ROW_FORMAT
/* Macros for getting the Octave(oct) ij that correspond to (row,col) [no pad involved] */
/* This one operates on GMT_MATRIX */
#define MEXM_IJ(M,row,col) ((row)*M->n_columns + (col))
/* And this on GMT_GRID */
#define MEXG_IJ(M,row,col) ((row)*M->header->nx + (col))
#else	/* Here we go for Matlab or Octave(mex) */
#ifdef GMT_MATLAB
#define MEX_PROG "Matlab"
#else
#define MEX_PROG "Octave(oct)"
#endif
#define MEX_COL_ORDER GMT_IS_COL_FORMAT
/* Macros for getting the Matlab/Octave(mex) ij that correspond to (row,col) [no pad involved] */
/* This one operates on GMT_MATRIX */
#define MEXM_IJ(M,row,col) ((col)*M->n_rows + (row))
/* And this on GMT_GRID */
#define MEXG_IJ(M,row,col) ((col)*M->header->ny + M->header->ny - (row) - 1)
#endif

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#define MODULE_LEN 32	/* Max length of a module name */

struct GMTMEX {	/* Array to hold information relating to output from GMT */
	enum GMT_enum_family type;	/* type of GMT data, i.e., GMT_IS_DATASET, GMT_IS_GRID, etc. */
	enum GMT_enum_std direction;	/* Either GMT_IN or GMT_OUT */
	int ID;				/* Registration ID returned by GMT_Register_IO */
	int lhs_index;			/* Corresponding index into plhs array */
	void *obj;			/* Pointer to this object */
};

EXTERN_MSC int GMTMEX_print_func (FILE *fp, const char *message);
EXTERN_MSC int GMTMEX_pre_process (void *API, const char *module, mxArray *plhs[], int nlhs, const mxArray *prhs[], int nrhs, const char *keys, struct GMT_OPTION **head, struct GMTMEX **X);
EXTERN_MSC int GMTMEX_post_process (void *API, struct GMTMEX *X, int n_items, mxArray *plhs[]);

#endif
