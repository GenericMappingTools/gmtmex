/*--------------------------------------------------------------------
 *	$Id$
 *
 *	Copyright (c) 1991-2015 by P. Wessel, W. H. F. Smith, R. Scharroo, J. Luis and F. Wobbe
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
 * 1) Create a new session and optionally return the API pointer. We provide for
 *    storing the pointer as a global variable (persistent) bewteen calls.
 * 2) Destroy a GMT session, either given the API pointer or by fetching it from
 *    the global (persistent) variable.
 * 3) Call any of the GMT modules.
 *
 * First argument to the gmt app is the API pointer, but it is optional once created.
 * Next argument is the command string that starts with the module name
 * Finally, there are optional comma-separated Matlab array entities required by the command.
 * Information about the options of each program is provided via GMT_Get_Moduleinfo.
 *
 * Version:	5.2
 * Created:	14-FEB-2015
 *
 */

#include "gmtmex.h"

#define MODULE_LEN 32

/* Being declared external we can access it between MEX calls */
static uintptr_t *pPersistent;    /* To store API address back and forth to a Matlab session */

/* Here is the exit function, which gets run when the MEX-file is
   cleared and when the user exits MATLAB. The mexAtExit function
   should always be declared as static. */
static void force_Destroy_Session(void) {
	void *API = (void *)pPersistent[0];	/* Get the GMT API pointer */
	if (API != NULL) {		/* Otherwise just silently ignore this call */
		mexPrintf("Destroying session due to a brute user usage.\n");
		if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT5 session\n");
	}
}

void usage(int nlhs, int nrhs) {
	/* Basic usage message */
	if (nrhs == 0) {	/* No arguments at all results in the GMT banner message */
		mexPrintf("\nGMT - The Generic Mapping Tools, Version %s\n", "5.2");
		mexPrintf("Copyright 1991-2015 Paul Wessel, Walter H. F. Smith, R. Scharroo, J. Luis, and F. Wobbe\n\n");
		mexPrintf("This program comes with NO WARRANTY, to the extent permitted by law.\n");
		mexPrintf("You may redistribute copies of this program under the terms of the\n");
		mexPrintf("GNU Lesser General Public License.\n");
		mexPrintf("For more information about these matters, see the file named LICENSE.TXT.\n");
		mexPrintf("For a brief description of GMT modules, type gmt ('--help')\n\n");
	}
	else {
		mexPrintf("Usage is:\n\tgmt ('create');  %% Create a new GMT/MEX session\n");
		mexPrintf("\tgmt ('module_name and options'[, <matlab arrays>]); %% Run a GMT module\n");
		mexPrintf("\tgmt ('destroy');  %% Destroy the GMT/MEX session\n");
		if (nlhs != 0)
			mexErrMsgTxt ("But meanwhile you already made an error by asking help and an output.\n");
	}
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int status = 0;                 /* Status code from GMT API */
	unsigned int first = 0;         /* Array ID of first command argument (not 0 when API-ID is first) */
	unsigned int help = 0;          /* 1 if we just gave --help */
	unsigned int got_API_in_input = 0; /* It will be set to 1 by gmt(API, 'module ...'); */
	int n_items = 0, pos = 0;       /* Number of Matlab arguments (left and right) */
	size_t str_length = 0, k = 0;   /* Misc. counters */
	struct GMTAPI_CTRL *API = NULL;	/* GMT API control structure */
	struct GMT_OPTION *options = NULL; /* Linked list of module options */
	struct GMTMEX *X = NULL;        /* Array of information about Matlab args */
	char *cmd = NULL;               /* Pointer used to get the user's Matlab command */
	char *opt_args = NULL;		/* Pointer to the user's module options */
	const char *keys = NULL;	/* This module's option keys */
	char module[MODULE_LEN];        /* Name of GMT module to call */
	uintptr_t *pti = NULL;          /* To locally store the API address */

	if (nrhs == 0) {	/* No arguments at all results in the GMT banner message */
		usage(nlhs, nrhs);
		return;
	}

	/* First check for the special commands create or destroy, while watching out for the lone --help argument */
	
	if (nrhs == 1) {	/* This may be create or --help */
		cmd = mxArrayToString (prhs[0]);
		help = !strncmp (cmd, "--help", 6U);
		if (help) {
			usage(nlhs, 1);
			return;
		}

		if (!strncmp (cmd, "create", 6U)) {
			if (pPersistent)                        /* See if have an GMT API pointer */
				API = (void *)pPersistent[0];
			if (API != NULL) {                      /* If another session still exists */
				mexPrintf ("A previous session is still active. Ignoring this 'create' request.\n");
				if (nlhs) {
					plhs[0] = mxCreateNumericMatrix (1, 0, mxUINT64_CLASS, mxREAL);
				}
				return;
			}

			/* Initializing new GMT session with zero pad and replacement printf function */
			if ((API = GMT_Create_Session ("GMT5", 0U, 3U, GMTMEX_print_func)) == NULL)
				mexErrMsgTxt ("Failure to create GMT5 Session\n");

			pPersistent = mxMalloc(sizeof(uintptr_t));
			pPersistent[0] = (uintptr_t)(API);
			mexMakeMemoryPersistent(pPersistent);

			if (nlhs) {	/* Return the API adress as an integer */
				plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
				pti = mxGetData(plhs[0]);
				*pti = *pPersistent;
			}

			mexAtExit(force_Destroy_Session);	/* Register an exit function. */
			return;
		}

		/* OK, no create and no --help, so it must be a single command with no arguments, nor the API. So get it */
		if (!pPersistent)
			mexErrMsgTxt ("Booo: you shouldn't have cleared this mex. Now the GMT5 session is lost (mem leaked).\n"); 
		API = (void *)pPersistent[0];	/* Get the GMT API pointer */
		if (API == NULL) mexErrMsgTxt ("This GMT5 session has already been destroyed, or currupted.\n"); 
		 
	}
	else if (mxIsScalar_(prhs[0]) && mxIsUint64(prhs[0])) {
		/* Here, nrhs > 1 . If first arg is a scalar int, assume it is the API memory adress */
		pti = (uintptr_t *)mxGetData(prhs[0]);
		API = (void *)pti[0];	/* Get the GMT API pointer */
		first = 1;		/* Commandline args start at prhs[1]. prhs[0] has the API id argument */
		got_API_in_input = 1;
	}
	else {		/* We still don't have the API */
		if (!pPersistent)
			mexErrMsgTxt("Booo: you shouldn't have cleared this mex. Now the GMT5 session is lost (mem leaked).\n"); 
		API = (void *)pPersistent[0];			/* Get the GMT API pointer */
		if (API == NULL) mexErrMsgTxt ("This GMT5 session has already been destroyed, or currupted.\n"); 
	}

	if (!cmd) 	/* First argument is the command string, e.g., 'blockmean -R0/5/0/5 -I1 or just destroy|free' */
		cmd = mxArrayToString (prhs[first]);

	/* WE CAN ALSO DESTROY THE SESSION BY SIMPLY CALLING "gmt('destroy')" */
	if (!strncmp (cmd, "destroy", 7U)) {
		if (nlhs != 0)
			mexErrMsgTxt ("Usage is gmt ('destroy');\n");

		if (GMT_Destroy_Options (API, &options)) mexErrMsgTxt ("Failure to destroy GMT5 options\n");
		if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT5 session\n");
		*pPersistent = 0;
		return;
	}

	/* Here we have a GMT module call of various sorts */
	
	/* 2. Get mex arguments, if any, and extract the GMT module name */
	str_length = strlen (cmd);				/* Length of command argument */
	for (k = 0; k < str_length && cmd[k] != ' '; k++);	/* Determine first space in command */
	memset (module, 0, MODULE_LEN * sizeof (char));		/* Initialize module name to blank */
	strncpy (module, cmd, k);				/* Isolate the module name in this string */

	/* 3. Determine the GMT module option keys, or list module usages and return if the module is not found.
	 *    Note: If module needed "gmt" prepended then module is returned with the correct prefix */
	if ((keys = GMT_Get_Moduleinfo (API, module)) == NULL) {
		GMT_Call_Module (API, NULL, GMT_MODULE_PURPOSE, NULL);
		if (GMT_Destroy_Session (API)) mexErrMsgTxt("Failure to destroy GMT5 session\n");
	}

	/* 4. Convert mex command line arguments to a linked option list */
	while (cmd[k] == ' ') k++;	/* Skip any spaces between module name and start of options */
	opt_args = (cmd[k]) ? &cmd[k] : NULL;
	if (opt_args && (options = GMT_Create_Options (API, 0, opt_args)) == NULL)
		mexErrMsgTxt ("Failure to parse GMT5 command options\n");

	/* 5. Parse the mex command, update GMT option lists, and register in/out resources, and return X array */
	pos = (got_API_in_input) ? 2 : 1;	/* Start depends on whether or not the API pointer was given or not */
	if ((n_items = GMTMEX_pre_process (API, module, plhs, nlhs, &prhs[MIN(pos,nrhs-1)], nrhs-pos, keys, &options, &X)) < 0)
		mexErrMsgTxt ("Failure to parse mex command options\n");
	
	/* 6. Run GMT module; give usage message if errors arise during parsing */
	status = GMT_Call_Module (API, module, GMT_MODULE_OPT, options);

	/* 7. Hook up module outputs to Matlab plhs arguments */
	if (GMTMEX_post_process (API, X, n_items, plhs)) mexErrMsgTxt ("Failure to extract GMT5-produced data\n");
	
	/* 8. Destroy linked option list */
	if (GMT_Destroy_Options (API, &options)) mexErrMsgTxt ("Failure to destroy GMT5 options\n");
}
