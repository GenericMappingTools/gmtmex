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

#include "gmtmex.h"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int status = 0;             /* Status code from GMT API */
	int module_id;              /* Module ID */
	int n_items = 0;            /* Number of Matlab arguments (left and right) */
	size_t str_length, k;       /* Misc. counters */
	struct GMTAPI_CTRL *API = NULL;     /* GMT API control structure */
	struct GMT_OPTION *options = NULL;  /* Linked list of options */
	struct GMTMEX *X = NULL;    /* Array of information about Matlab args */
	char *cmd = NULL;           /* Pointer used to get Matlab command */
	char module[BUFSIZ];        /* Name of GMT module to call */

	if (nrhs == 0) {
		mexPrintf("\ngmt - The Generic Mapping Tools, Version %s\n", "5.0");
		mexPrintf("Copyright 1991-2013 Paul Wessel, Walter H. F. Smith, R. Scharroo, J. Luis, and F. Wobbe\n\n");
		mexPrintf("This program comes with NO WARRANTY, to the extent permitted by law.\n");
		mexPrintf("You may redistribute copies of this program under the terms of the\n");
		mexPrintf("GNU Lesser General Public License.\n");
		mexPrintf("For more information about these matters, see the file named LICENSE.TXT.\n");
		mexPrintf("For a brief description of GMT programs, type gmt5('--help')\n\n");
		return;
	}

	/* 1. Initializing new GMT session with zero pad and replacement printf function */
	if ((API = GMT_Create_Session ("GMT5", 0U, 1U, GMTMEX_print_func)) == NULL) mexErrMsgTxt ("Failure to create GMT5 Session\n");

	/* 2. Get mex arguments, if any, and extract the GMT module name */
	cmd = mxArrayToString (prhs[0]);			/* First argument is the command string, e.g., 'blockmean -R0/5/0/5 -I1' */
	str_length = strlen  (cmd);				/* Length of command argument */
	for (k = 0; k < str_length && cmd[k] != ' '; k++);	/* Determine first space in command */
	memset ((void *)module, 0, BUFSIZ*sizeof (char));	/* Initialize module name to blank */
	strncpy (module, cmd, k);				/* Isolate the module name in this string */

	/* 3. Determine the GMT module ID, or list module usages and return if module is not found */
	if ((module_id = gmtmex_find_module (module)) == -1) {
		GMT_List_Module (API, NULL);
		if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT5 session\n");
		return;
	}

	/* 4. Convert mex command line arguments to a linked option list */
	if ((options = GMT_Create_Options (API, 0, &cmd[k+1])) == NULL)
		mexErrMsgTxt ("Failure to parse GMT5 command options\n");

	/* 5. Parse the mex command, update GMT option lists, and register in/out resources, and return X array */
	if ((n_items = GMTMEX_pre_process (API, module, plhs, nlhs, &prhs[1], nrhs-1, keys[module_id], options, &X)) < 0)
		mexErrMsgTxt ("Failure to parse mex command options\n");
	
	/* 6. Run GMT module; give usage message if errors arise during parsing */
	status = GMT_Call_Module (API, module, -1, options);

	/* 7. Hook up module output to Matlab plhs arguments */
	if (GMTMEX_post_process (API, X, n_items, plhs)) mexErrMsgTxt ("Failure to extract GMT5-produced data\n");
	
	/* 8. Destroy linked option list */
	if (GMT_Destroy_Options (API, &options)) mexErrMsgTxt ("Failure to destroy GMT5 options\n");

	/* 9. Destroy GMT session */
	if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT5 session\n");
}
