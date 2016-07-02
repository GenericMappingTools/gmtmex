/*
 *	$Id$
 *
 *	Copyright (c) 1991-2016 by P. Wessel and J. Luis
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
/* GMT convenience functions used by MATLAB/OCTAVE mex/oct API
 */

#ifndef GMTMEX_H
#define GMTMEX_H

#include "gmt.h"
#include <string.h>
#include <ctype.h>
#include <limits.h>

#ifdef GMT_OCTOCT
#	include <oct.h>
#else
#	include <mex.h>
#	define mxIsScalar_(mx) \
		( (2 == mxGetNumberOfDimensions(mx)) \
		&&  (1 == mxGetM(mx))&&  (1 == mxGetN(mx)) )
#endif	/* Matlab and Octave(mex) */

/* Matlab and Octave (in -mex mode) are identical, oct files are different and not yet tested */

/* Older Ml versions don't have mwSize */
#ifndef GMT_OCTMEX
#if !defined(MATLAB_VERSION)
#if !defined(MWSIZE_MAX)
#    define MATLAB_VERSION 0x2006a /* R2006a or earlier */
#elif MX_API_VER < 0x07040000
#    define MATLAB_VERSION 0x2006b /* R2006b */
#elif !defined(FMT_PTRDIFF_T)
#    define MATLAB_VERSION 0x2007a /* R2007a */
#elif !defined(CUINT64_T)
#    define MATLAB_VERSION 0x2007b /* R2007b */
#elif defined(mxSetLogical)
#    define MATLAB_VERSION 0x2008a /* R2008a */
#else
#    if !defined(blas_h)
#        include "blas.h"
#    endif
#    if !defined(lapack_h)
#        include "lapack.h"
#    endif
#    if !defined(MATHWORKS_MATRIX_MATRIX_PUB_FWD_H)
#        if defined(CHAR16_T)
#            if !defined(COMPLEX_TYPES)
#                define MATLAB_VERSION 0x2008b /* R2008b */
#            elif !defined(cgeqr2p)
#                define MATLAB_VERSION 0x2010b /* R2010b */
#            else
#                define MATLAB_VERSION 0x2011a /* R2011a */
#            endif
#        else
#            include "emlrt.h"
#            define MATLAB_VERSION EMLRT_VERSION_INFO /* R2011b or later */
#        endif
#    else
#        if !defined(COMPLEX_TYPES)
#            define MATLAB_VERSION 0x2009a /* R2009a */
#        elif !defined(cgbequb)
#            define MATLAB_VERSION 0x2009b /* R2009b */
#        else
#            define MATLAB_VERSION 0x2010a /* R2010a */
#        endif
#    endif
#endif
#endif /* if !defined(MATLAB_VERSION) */

#if MATLAB_VERSION < 0x2006b
typedef int mwSize;
#endif
#endif

#ifdef GMT_OCTOCT	/* Octave oct files only */
#	define MEX_PROG "Octave(oct)"
#	define MEX_COL_ORDER GMT_IS_ROW_FORMAT
	/* Macros for getting the Octave(oct) ij that correspond to (row,col) [no pad involved] */
	/* This one operates on GMT_MATRIX */
#	define MEXM_IJ(M,row,col) ((row)*M->n_columns + (col))
	/* And this on GMT_GRID */
#	define MEXG_IJ(M,row,col) ((row)*M->header->nx + (col))
#else	/* Here we go for Matlab or Octave(mex) */
#	ifdef GMT_MATLAB
#		define MEX_PROG "Matlab"
#	else
#		define MEX_PROG "Octave(mex)"
#	endif
#	define MEX_COL_ORDER GMT_IS_COL_FORMAT
	/* Macros for getting the Matlab/Octave(mex) ij that correspond to (row,col) [no pad involved] */
	/* This one operates on GMT_MATRIX */
#	define MEXM_IJ(M,row,col) ((col)*M->n_rows + (row))
	/* And this on GMT_GRID */
#	define MEXG_IJ(M,row,col) ((col)*M->header->ny + M->header->ny - (row) - 1)
#endif

/* Macro for indecing into a GMT grid [with pad] */
#define GMT_IJP(h,row,col) ((uint64_t)(((int64_t)(row)+(int64_t)h->pad[GMT_YHI])*((int64_t)h->mx)+(int64_t)(col)+(int64_t)h->pad[GMT_XLO]))

#define MODULE_LEN 	32	/* Max length of a GMT module name */

EXTERN_MSC int GMTMEX_print_func        (FILE *fp, const char *message);
EXTERN_MSC void * GMTMEX_Get_Grid       (void *API, struct GMT_GRID *G);
EXTERN_MSC void * GMTMEX_Get_Dataset    (void *API, struct GMT_DATASET *D);
EXTERN_MSC void * GMTMEX_Get_Textset    (void *API, struct GMT_TEXTSET *M);
EXTERN_MSC void * GMTMEX_Get_CPT        (void *API, struct GMT_PALETTE *P);
EXTERN_MSC void * GMTMEX_Get_Image      (void *API, struct GMT_IMAGE *I);
EXTERN_MSC void * GMTMEX_Register_IO    (void *API, struct GMT_RESOURCE *X, const mxArray *ptr);
EXTERN_MSC void GMTMEX_Free_Textset     (void *API, struct GMT_TEXTSET *T);
#if GMT_MINOR_VERSION > 2
EXTERN_MSC void * GMTMEX_Get_POSTSCRIPT (void *API, struct GMT_POSTSCRIPT *P);
#endif
mxClassID GMTMEX_type (void *API);
#endif
#if 0
EXTERN_MSC void * GMTMEX_Get_Dataset    (void *API, struct GMT_VECTOR *V);
#endif
