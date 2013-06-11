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
 * This is the Matlab/Octave GMT application, which can do the following:
 * 1) Create a new session and return the API pointer.
 * 2) Destroy a GMT session given the API pointer
 * 3) Call any of the GMT modules.
 * First argument to GMT is expected to be the API, followed by a command
 * string, with optional comma-separated Matlab array entities.
 * Information about the options of each program is provided by the include
 * files generated from mexproginfo.txt.
 *
 * Version:	5
 * Created:	12-May-2013
 *
 */

#include "gmtmex.h"

/* Being declared external we can access it between MEX calls */
static uintptr_t *pti, *pPersistent;    /* To store API address back and forth to a Matlab session */

/* Here is the exit function, which gets run when the MEX-file is
   cleared and when the user exits MATLAB. The mexAtExit function
   should always be declared as static. */
static void force_Destroy_Session(void) {
	void *API;
	mexPrintf("Destroying session due to a brute user usage.\n");
	API = (void *)pPersistent[0];			/* Get the GMT API pointer */
	if (API == NULL) mexErrMsgTxt ("Grrr: this GMT5 session has already been destroyed.\n"); 
	if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT5 session\n");
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int status = 0;                 /* Status code from GMT API */
	unsigned int first;             /* Array ID of first command argument */
	bool help;                      /* True if we just gave --help */
	int n_items = 0;                /* Number of Matlab arguments (left and right) */
	int module_id;
	size_t str_length, k;           /* Misc. counters */
	struct GMTAPI_CTRL *API = NULL;	/* GMT API control structure */
	struct GMT_OPTION *options = NULL; /* Linked list of options */
	struct GMTMEX *X = NULL;        /* Array of information about Matlab args */
	char *cmd = NULL;               /* Pointer used to get Matlab command */
	char module[BUFSIZ];            /* Name of GMT module to call */

	if (nrhs == 0) {	/* No arguments at all results in the GMT banner message */
		mexPrintf("\nGMT - The Generic Mapping Tools, Version %s\n", "5.0");
		mexPrintf("Copyright 1991-2013 Paul Wessel, Walter H. F. Smith, R. Scharroo, J. Luis, and F. Wobbe\n\n");
		mexPrintf("This program comes with NO WARRANTY, to the extent permitted by law.\n");
		mexPrintf("You may redistribute copies of this program under the terms of the\n");
		mexPrintf("GNU Lesser General Public License.\n");
		mexPrintf("For more information about these matters, see the file named LICENSE.TXT.\n");
		mexPrintf("For a brief description of GMT modules, type GMT('--help')\n\n");
		return;
	}

	/* First check for the special commands create or destroy, while watching out for the lone --help argument */
	
	if (nrhs == 1) {	/* This may be create or --help */
		cmd = mxArrayToString (prhs[0]);
		help = !strncmp (cmd, "--help", 6U);
		if (!strncmp (cmd, "create", 6U) || help) {
			if (!help && nlhs != 1) {
				mexErrMsgTxt ("Usage is API = GMT ('create');\n");
			}
			/* Initializing new GMT session with zero pad and replacement printf function */
			if ((API = GMT_Create_Session ("GMT5", 0U, 1U, GMTMEX_print_func)) == NULL) mexErrMsgTxt ("Failure to create GMT5 Session\n");
			if (!help) {	/* This was create, so just return API pointer */
				plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
				pti = mxGetData(plhs[0]);
				pti[0] = (uintptr_t)(API);
				pPersistent = mxMalloc(sizeof(uintptr_t));
				pPersistent = pti;
				mexMakeMemoryPersistent(pPersistent);
				mexAtExit(force_Destroy_Session);	/* Register an exit function. */
				return;
			}
			first = 0;		/* Commandline args start at prhs[0] since there is no API argument */
		}
	}
	else {	/* Here, nrhs > 1 */
		pti = (uintptr_t *)mxGetData(prhs[0]);
		API = (void *)pti[0];	/* Get the GMT API pointer */
		first = 1;		/* Commandline args start at prhs[1] */
	}

	/* WE CAN ALSO DESTROY THE SESSION BY SIMPLY CALLING "gmt('destroy')" */
	if (nlhs == 0 && nrhs == 1 && !strncmp (cmd, "destroy", 7U)) {	/* Destroy GMT session */
		if (!pti) mexErrMsgTxt ("Booo: you shouldn't have cleared this mex. Now the GMT5 session is lost (mem leaked).\n"); 
		API = (void *)pti[0];			/* Get the GMT API pointer */
		if (API == NULL) mexErrMsgTxt ("Grrr: this GMT5 session has already been destroyed.\n"); 
		if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT5 session\n");
		*pti = 0;
		return;
	}

	cmd = mxArrayToString (prhs[first]);	/* First argument is the command string, e.g., 'blockmean -R0/5/0/5 -I1 or just destroy|free' */
	if (!strncmp (cmd, "destroy", 7U)) {	/* Destroy GMT session */
		if (!(nlhs == 0 && nrhs == 2)) {
			mexErrMsgTxt ("Usage is GMT (API, 'destroy');\n");
		}
		if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT5 session\n");
		return;
	}

	/* Here we have GMT module calls of various sorts */
	
	/* 2. Get mex arguments, if any, and extract the GMT module name */
	str_length = strlen  (cmd);				/* Length of command argument */
	for (k = 0; k < str_length && cmd[k] != ' '; k++);	/* Determine first space in command */
	memset ((void *)module, 0, BUFSIZ*sizeof (char));	/* Initialize module name to blank */
	strncpy (module, cmd, k);				/* Isolate the module name in this string */

	/* 3. Determine the GMT module ID, or list module usages and return if module is not found */
	if ((module_id = gmtmex_find_module (module)) == -1) {
		GMT_List_Module (API, NULL);
		return;
	}

	/* 4. Convert mex command line arguments to a linked option list */
	if ((options = GMT_Create_Options (API, 0, &cmd[k+1])) == NULL)
		mexErrMsgTxt ("Failure to parse GMT5 command options\n");

	/* 5. Parse the mex command, update GMT option lists, and register in/out resources, and return X array */
	if ((n_items = GMTMEX_pre_process (API, module, plhs, nlhs, &prhs[2], nrhs-2, keys[module_id], options, &X)) < 0)
		mexErrMsgTxt ("Failure to parse mex command options\n");
	
	/* 6. Run GMT module; give usage message if errors arise during parsing */
	status = GMT_Call_Module (API, module, -1, options);

	/* 7. Hook up module output to Matlab plhs arguments */
	if (GMTMEX_post_process (API, X, n_items, plhs)) mexErrMsgTxt ("Failure to extract GMT5-produced data\n");
	
	/* 8. Destroy linked option list */
	if (GMT_Destroy_Options (API, &options)) mexErrMsgTxt ("Failure to destroy GMT5 options\n");
}

