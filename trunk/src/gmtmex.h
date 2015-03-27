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

#ifdef GMT_OCTOCT
#include <oct.h>
#else
#include <mex.h>
#define mxIsScalar_(mx) \
	( (2 == mxGetNumberOfDimensions(mx)) \
	&&  (1 == mxGetM(mx))&&  (1 == mxGetN(mx)) )
#endif	/* Matlab and Octave(mex) */

/* Matlab and Octave (in -mex mode) are identical, oct files are different and not yet tested */

/* Older Ml versions don't have mwSize */
#ifndef GMT_OCTMEX
#ifndef mwSize
	typedef int mwSize;
#endif
#endif

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
#define MEX_PROG "Octave(mex)"
#endif
#define MEX_COL_ORDER GMT_IS_COL_FORMAT
/* Macros for getting the Matlab/Octave(mex) ij that correspond to (row,col) [no pad involved] */
/* This one operates on GMT_MATRIX */
#define MEXM_IJ(M,row,col) ((col)*M->n_rows + (row))
/* And this on GMT_GRID */
#define MEXG_IJ(M,row,col) ((col)*M->header->ny + M->header->ny - (row) - 1)
#endif

#define MODULE_LEN 	32	/* Max length of a GMT module name */
#define ARG_MARKER	'$'	/* Character that indicates an memory reference to data */

EXTERN_MSC int GMTMEX_print_func (FILE *fp, const char *message);
EXTERN_MSC void * GMTMEX_Get_Grid    (void *API, struct GMT_GRID *G);
EXTERN_MSC void * GMTMEX_Get_Dataset (void *API, struct GMT_MATRIX *M);
EXTERN_MSC void * GMTMEX_Get_Textset (void *API, struct GMT_TEXTSET *M);
EXTERN_MSC void * GMTMEX_Get_CPT     (void *API, struct GMT_PALETTE *P);
EXTERN_MSC void * GMTMEX_Get_Image   (void *API, struct GMT_IMAGE *I);
EXTERN_MSC void * GMTMEX_Register_IO (void *API, unsigned int data_type, unsigned int geometry, unsigned int direction, const mxArray *ptr, int *ID);
EXTERN_MSC void GMTMEX_Free_Textset (void *API, struct GMT_TEXTSET *T);
#endif
