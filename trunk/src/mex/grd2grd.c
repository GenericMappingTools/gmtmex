/*
 *	$Id: grd2grd.c,v 1.2 2011-03-15 02:06:37 guru Exp $
 *
 *	Copyright (c) 1991-2011 by P. Wessel, W. H. F. Smith, R. Scharroo, and J. Luis
 *      See LICENSE.TXT file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 of the License.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: www.soest.hawaii.edu/pwessel
 *--------------------------------------------------------------------*/
/* Program:	MATLAB/OCTAVE interface to modules that take one grid in
 *			      and create another grid in output
 */
 
#include "gmt_mex.h"

/* EXTERN_MSC GMT_LONG FUNC (struct GMTAPI_CTRL *API, GMT_LONG argc, char **cmd); */

/* Matlab Gateway routine for this grd??? prog */

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	GMT_LONG status;
	struct GMTAPI_CTRL *API = NULL;		/* GMT API control structure */
	struct GMT_GRID *Gin = NULL, *Gout = NULL;
	char *input = NULL, *output = NULL, *options = NULL, *cmd = NULL; 

	/* Make sure in/out arguments match expectation, or give usage message */
	if ((nrhs < 2 || nrhs > 5) || nlhs > 2) {	/* Not what we expected, or nothing */
		GMT5MEX_banner;
		mexPrintf ("usage: Z = grdXXX ('filename', 'options');\n");
		mexPrintf ("	Z = grdXXX (Z, info, 'options');\n");
		mexPrintf ("	[Z info] = grdXXX ('filename', 'options');\n");
		mexPrintf ("	[Z info] = grdXXX (Z, info, 'options');\n");
		return;
	}
	if (!mxIsChar(prhs[nrhs-1])) mexErrMsgTxt ("Last input must contain the options string\n");

	/* Initializing new GMT session */
	if (GMT_Create_Session (&API, "MEX", FUNC_MODE)) mexErrMsgTxt ("Failure to create GMT Session\n");

	/* Make sure options are given, and get them */
	options = GMTMEX_options_init (API, prhs, nrhs);

	/* Set up input grid (actual or via Matlab matriz) */
	input = GMTMEX_src_grid_init (API, prhs, nrhs, &Gin);

	/* Register a grid struct G to be the destination, allocated and written to by the module */
	output = GMTMEX_dest_grid_init (API, &Gout, nlhs, options);

	/* Build module command from input, ouptput, and option strings */
	cmd = GMTMEX_build_cmd (API, input, options, output, GMT_IS_GRID);

	/* Run GMT_grdXXX module, or give usage message if errors arise during parsing */
	if ((status = FUNC (API, 0, (void *)cmd))) mexErrMsgTxt ("Run-time error\n");
	
	/* Pass output arguments to Matlab vectors Z, with optional (x, y) or hdr. */
	if (nlhs) GMTMEX_prep_mexgrd (API, plhs, nlhs, Gout);

	/* Destroy the main grid returned from the module  */
	GMT_Destroy_Data (API, GMT_REFERENCE, (void **)&Gout);
	
	/* Free temporary local variables  */
	GMTMEX_free (input, output, options, cmd);
	GMT_free_grid (API->GMT, &Gin, TRUE);	/* TRUE since we had to duplicate the Matlab matrix */
	
	/* Destroy GMT API session */
	if (GMT_Destroy_Session (&API)) mexErrMsgTxt ("Failure to destroy GMT Session\n");
}
