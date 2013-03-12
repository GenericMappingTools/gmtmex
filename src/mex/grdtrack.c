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
/* Program:	MATLAB/OCTAVE interface to grdtrack
 */
 
#include "gmt_mex.h"

/* Matlab Gateway routine for grdtrack */

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int status;
	struct	GMTAPI_CTRL *API = NULL;		/* GMT API control structure */
	struct	GMT_VECTOR *Vi = NULL, *Vo = NULL;
	float	*Z = NULL;
	char	*input = NULL, *inputG = NULL, *output = NULL, *options = NULL, *cmd = NULL; 
	int	n_cols, n_start;

	/* Make sure in/out arguments match expectation, or give usage message */
	if (!(nrhs == 1 || (nrhs >= 3 || nrhs <= 5)) || !(nlhs == 0 || nlhs == 1)) {
		GMT5MEX_banner;
		mexPrintf ("usage: z = grdtrack (Z, hdr, x, y, 'options');\n");
		mexPrintf ("	   z = grdtrack (Z, hdr, x, y);\n");
		mexPrintf ("	   z = grdtrack (x, y, 'options');\n");
		mexPrintf ("	   z = grdtrack (Z, hdr, 'options');\n");
		mexPrintf ("	   z = grdtrack ('options');\n");
		mexPrintf ("	       grdtrack ('options')\n");
		return;
	}
	if (!mxIsChar(prhs[nrhs-1]) && nrhs != 4) mexErrMsgTxt ("Last input must contain the options string\n");

	/* Initializing new GMT session */
	if ((API = GMT_Create_Session ("GMT/MEX-API", 2U, 0U)) == NULL) mexErrMsgTxt ("Failure to create GMT Session\n");

	/* Make sure options are given, and get them */
	options = GMTMEX_options_init (API, prhs, nrhs);

	/* Set up input grid (actual or via Matlab matrix) */
	if (!mxIsChar(prhs[0]) && mxGetM(prhs[0]) > 1 && mxGetN(prhs[0]) > 1) 
		inputG = GMTMEX_src_grid_init (API, prhs, nrhs);

	/* Set up input file (actual or via Matlab vectors) */
	n_cols = 2;	n_start = 2;
	if (nrhs == 3 && (mxGetM(prhs[0]) == 1 || mxGetM(prhs[0]) == 1))	/* grdtrack (x, y, 'options') */
		n_start = 0;
	else if (nrhs == 1 || nrhs == 3) {		/* No data vectors in input */
		n_cols = n_start = 0;
	}
	input = GMTMEX_src_vector_init (API, prhs, n_cols, n_start, &Vi);

	/* Cat the two input strings (we can only send one to GMTMEX_build_cmd()) */
	if (inputG) strcat(input, inputG);

	/* Register output vectors Vo to be the destination, allocated and written to by the module */
	output = GMTMEX_dest_vector_init (API, (int)nlhs, &Vo, nlhs, options);

	/* Build module command from input, ouptput, and option strings */
	cmd = GMTMEX_build_cmd (API, input, options, output, GMT_IS_DATASET);
	mexPrintf ("cmd = (%s)\n", cmd);

	/* Run blockmean module, and give usage message if errors arise during parsing */
	if ((status = GMT_grdtrack (API, 0, cmd))) mexErrMsgTxt ("Run-time error\n");
	
	/* Pass output arguments to Matlab column vectors. */
	if (nlhs) GMTMEX_prep_mextbl (API, plhs, nlhs, Vo);

	/* Destroy the columns returned from the module  */
	GMT_free_vector (API->GMT, &Vo, true);	/* true since vectors are being duplicated for Matlab */
	
	/* Free temporary local variables  */
	GMTMEX_free (input, output, options, cmd);
	GMT_free_vector (API->GMT, &Vi, false);	/* false since vectors came from Matlab */
	
	/* Destroy GMT API session */
	if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT Session\n");

	return;
}
