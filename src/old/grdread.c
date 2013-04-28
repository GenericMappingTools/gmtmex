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
/* Program:	grdread.c
 * Purpose:	matlab/octave callable routine to read a GMT grid file
 */
 
#include "gmt_mex.h"

/* Matlab Gateway routine */

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	struct GMTAPI_CTRL *API = NULL;		/* GMT API control structure */
	struct GMT_GRID *G = NULL;
	float	*z = NULL;
	char *filein = NULL;
	int row, col;
	uint64_t gmt_node;
	int	px = -1, py = -1, pz = -1, pi = -1;

	if (nrhs != 1 || nlhs < 1 || nlhs > 4) {	/* Give usage message and return */
		GMT5MEX_banner;
		mexPrintf ("usage: z = grdread ('filename');\n");
		mexPrintf ("	[z info] = grdread ('filename');\n");
		mexPrintf ("	[x y z] = grdread ('filename');\n");
		mexPrintf ("	[x y z info] = grdread ('filename');\n");
		return;
	}
	if (!mxIsChar(prhs[nrhs-1])) mexErrMsgTxt ("Input must contain the filename string\n");

	/* Load the file name into a char string */

	filein = mxArrayToString (prhs[0]);	/* Load the file name into a char string */

	/* 1. Initializing new GMT session */
	if ((API = GMT_Create_Session ("GMT/MEX-API", 2U, 0U)) == NULL) mexErrMsgTxt ("Failure to create GMT Session\n");

	/* 2. READING IN A GRID */
	if ((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, filein, NULL)) == NULL)
		mexErrMsgTxt ("GMT: (grdread) Read failure\n");
	
	/* Create a matrix for the return array */

	pz = (nlhs >= 3) ? 2 : 0;
	if (nlhs == 2 || nlhs == 4) pi = nlhs - 1;
	if (nlhs > 2) {px = 0; py = 1;}

	plhs[pz] = mxCreateNumericMatrix (G->header->ny, G->header->nx, mxSINGLE_CLASS, mxREAL);
	z = mxGetData (plhs[pz]);
	
	/*  Load the real grd array into a double matlab array
	    by transposing from padded GMT grd format to unpadded matlab format */
    
	GMT_grd_loop (API->GMT, G, row, col, gmt_node) z[MEX_IJ(G,row,col)] = G->data[gmt_node];
	    
	/* Create scalars for return arguments */

	if (pi >= 0) GMTMEX_grdheader2info (plhs, G, pi);	/* Also return info array */
	if (px >= 0) GMTMEX_grdxy (API, plhs, G, px, py);	/* Return x,y arrays also */
	
	if (GMT_Destroy_Data (API, GMT_ALLOCATED, &G)) mexErrMsgTxt ("Run-time error\n");
	
	/* 10. Destroy GMT API session */
	if (GMT_Destroy_Session (API)) mexErrMsgTxt ("GMT: (surface) Failure to destroy GMT Session\n");

	return;
}