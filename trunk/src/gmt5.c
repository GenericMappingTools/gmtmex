/*--------------------------------------------------------------------
 *	$Id$
 *
 *	Copyright (c) 1991-$year by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis and F. Wobbe
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation; version 3 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * This is the Matlab/Octave gmt5 application, which can call any of the
 * GMT modules. First argument to gmt5 is expected to be the module name.
 * Information about the options of each program is provided by the include
 * files generated from mexproginfo.txt.
 *
 * Version:	5
 * Created:	20-Apr-2013
 *
 */

#include "gmt_mex.h"
#include "gmtmex_id.h"
#include "gmtmex_keys.h"
#include "gmtmex_progs.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int status = 0;				/* Status code from GMT API */
	int module_id;    /* Module ID */
	size_t str_length, k;
	struct GMTAPI_CTRL *API = NULL;		/* GMT API control structure */
	struct GMT_OPTION *options = NULL;	/* Linked list of options */
	char *cmd = NULL;
	char module[BUFSIZ];

	/* 1. Initializing new GMT session */
	if ((API = GMT_Create_Session ("GMT5", 2U, 1U, GMTMEX_print_func)) == NULL) mexErrMsgTxt ("Failure to create GMT Session\n");

	/* 2. Get mex arguments, if any, and extract the module name */
	cmd = mxArrayToString (prhs[0]);		/* First argument is the command string, e.g., 'blockmean $ -R0/5/0/5 -I1' */
	str_length = strlen  (cmd);			/* Length of command argument */
	for (k = 0; k < str_length && cmd[k] != ' '; k++);	/* Determine first space in command */
	memset ((void *)module, 0, BUFSIZ*sizeof (char));
	strncpy (module, cmd, k);			/* Isolate the module name */

	/* 3. Determine the GMT module ID, or list module usages and return if not found */
	if ((module_id = GMT_Get_Module (API, module)) == GMT_ID_NONE) {
		GMT_List_Module (API, -1);
		return;
	}

	/* 4. Convert mex command line arguments to a linked option list */
	if ((options = GMT_Create_Options (API, 0, &cmd[k+1])) == NULL) mexErrMsgTxt ("Failure to parse GMT command options\n");

	/* 5. Parse the mex command, update GMT option lists, register in/out resources */
	if (GMTMEX_parser (API, plhs, nlhs, prhs, nrhs, keys[module_id], options)) mexErrMsgTxt ("Failure to parse mex command options\n");
	
	/* 6. Run GMT cmd module, or give usage message if errors arise during parsing */
	status = GMT_Call_Module (API, module_id, -1, options);

	/* 7. Destroy linked option list */
	if (GMT_Destroy_Options (API, &options)) mexErrMsgTxt ("Failure to destroy GMT options\n");

	/* 8. Destroy GMT session */
	if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT session\n");
}
