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

#ifndef GMTMEX_H
#define GMTMEX_H

#include "gmt.h"
#include <mex.h>
#include <string.h>
#include <ctype.h>

#ifdef GMT_MATLAB
#define MEX_PROG "Matlab"
#else
#define MEX_PROG "Octave"
#endif

struct GMTMEX {	/* Array to hold information relating to output from GMT */
	unsigned int type;	/* type of GMT data, i.e., GMT_IS_DATASET, GMT_IS_GRID, etc. */
	int ID;			/* Registration ID returned by GMT_Register_IO */
	int lhs_index;		/* Corresponding index into plhs array */
};

extern int GMTMEX_print_func (FILE *fp, const char *message);
extern int GMTMEX_pre_process (void *API, mxArray *plhs[], int nlhs, const mxArray *prhs[], int nrhs, char *keys, struct GMT_OPTION *head, struct GMTMEX **X);
extern int GMTMEX_post_process (void *API, struct GMTMEX *X, int n_items, mxArray *plhs[]);

#endif
