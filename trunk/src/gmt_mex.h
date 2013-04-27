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

#ifndef GMT_MEX_H
#define GMT_MEX_H

#include "gmt.h"
#include <mex.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#if defined(WIN32) && !defined(lrint)
#	define lrint (int64_t)rint
#endif

#define GMT_VIA_MEX	0	/* See what this is for later */
#ifdef GMT_MATLAB
#define MEX_PROG "Matlab"
#else
#define MEX_PROG "Octave"
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

/* Macros for getting the Matlab/Octave ij that correspond to (row,col) [no pad involved] */
/* This one operates on GMT_MATRIX */
#define MEXM_IJ(M,row,col) ((col)*M->n_rows + M->n_rows - (row) - 1)

/* And this on GMT_GRID */
#define MEXG_IJ(M,row,col) ((col)*M->header->ny + M->header->ny - (row) - 1)

#define GMT5MEX_banner mexPrintf("The Generic Mapping Tools v. 5 %s interface\n", MEX_PROG)
#define GMT_IS_PS	99	/* Use for PS output; use GMT_IS_GRID or GMT_IS_DATASET for data */

struct GMTMEX {	/* Array to hold information relating to output from GMT */
	unsigned int type;	/* type of GMT data, i.e., GMT_IS_DATASET, GMT_IS_GRID, etc. */
	int ID;			/* Registration ID returned by GMT_Register_IO */
	int lhs_index;		/* Corresponding index into plhs array */
};

int GMTMEX_print_func (FILE *fp, const char *message);

void GMTMEX_grdheader2info (mxArray *plhs[], struct GMT_GRID *G, int item);
void GMTMEX_grdxy (void *API, mxArray *plhs[], struct GMT_GRID *G, int px, int py);
void GMTMEX_prep_mexgrd (void *API, mxArray *plhs[], int nlhs, struct GMT_GRID *G);
void GMTMEX_prep_mextbl (void *API, mxArray *plhs[], int nlhs, struct GMT_VECTOR *V);
double *GMTMEX_info2grdheader (void *API, const mxArray *prhs[], int nrhs, struct GMT_GRID_HEADER *header);
char *GMTMEX_src_grid_init (void *API, const mxArray *prhs[], int nrhs);
char *GMTMEX_src_vector_init (void *API, const mxArray *prhs[], unsigned int n_cols, int n_start, struct GMT_VECTOR **V);
char *GMTMEX_dest_grid_init (void *API, int *ID, int nlhs, char *options);
char *GMTMEX_dest_vector_init (void *API, unsigned int n_cols, struct GMT_VECTOR **V, int nlhs, char *options);
char *GMTMEX_options_init (void *API, const mxArray *prhs[], int nrhs);
char *GMTMEX_build_cmd (void *API, char *src, char *options, char *dest, int mode);
void GMTMEX_free (char *input, char *output, char *options, char *cmd);
int GMTMEX_pre_process (void *API, void *plhs[], int nlhs, void *prhs[], int nrhs, char *keys, struct GMT_OPTION *head, struct GMTMEX **X);
int GMTMEX_post_process (void *API, struct GMTMEX *X, int n_items, mxArray *plhs[]);

#endif
