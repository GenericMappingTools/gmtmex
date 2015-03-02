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
 * This is the Matlab/Octave(mex) GMT application, which can do the following:
 * 1) Create a new session and optionally return the API pointer. We provide for
 *    storing the pointer as a global variable (persistent) between calls.
 * 2) Destroy a GMT session, either given the API pointer or by fetching it from
 *    the global (persistent) variable.
 * 3) Call any of the GMT modules while passing data in and out of GMT.
 *
 * First argument to the gmt function is the API pointer, but it is optional once created.
 * Next argument is the command string that starts with the module name
 * Finally, there are optional comma-separated Matlab array entities required by the command.
 * Information about the options of each program is provided via GMT_Encode_Options.
 *
 * Version:	5.2
 * Created:	20-FEB-2015
 *
 */

#include "gmtmex.h"

extern int GMT_get_V (char arg);

/* Being declared external we can access it between MEX calls */
static uintptr_t *pPersistent;    /* To store API address back and forth within a single Matlab session */

/* Here is the exit function, which gets run when the MEX-file is
   cleared and when the user exits MATLAB. The mexAtExit function
   should always be declared as static. */
static void force_Destroy_Session (void) {
	void *API = (void *)pPersistent[0];	/* Get the GMT API pointer */
	if (API != NULL) {		/* Otherwise just silently ignore this call */
		mexPrintf("GMT: Destroying GMT session due to a brute user usage.\n");
		if (GMT_Destroy_Session (API)) mexErrMsgTxt ("Failure to destroy GMT session\n");
		mexPrintf("GMT: GMT session destroyed.\n");
	}
}

void usage (int nlhs, int nrhs) {
	/* Basic usage message */
	if (nrhs == 0) {	/* No arguments at all results in the GMT banner message */
		mexPrintf("\nGMT - The Generic Mapping Tools, Version %s %s API\n", "5.2", MEX_PROG);
		mexPrintf("Copyright 1991-2015 Paul Wessel, Walter H. F. Smith, R. Scharroo, J. Luis, and F. Wobbe\n\n");
		mexPrintf("This program comes with NO WARRANTY, to the extent permitted by law.\n");
		mexPrintf("You may redistribute copies of this program under the terms of the\n");
		mexPrintf("GNU Lesser General Public License.\n");
		mexPrintf("For more information about these matters, see the file named LICENSE.TXT.\n");
		mexPrintf("For a brief description of GMT modules, type gmt ('help')\n\n");
	}
	else {
		mexPrintf("Usage is:\n\tgmt ('create');  %% Create a new GMT/MEX session\n");
		mexPrintf("\tgmt ('module_name and options'[, <matlab arrays>]); %% Run a GMT module\n");
		mexPrintf("\tgmt ('destroy');  %% Destroy the GMT/MEX session\n");
		if (nlhs != 0)
			mexErrMsgTxt ("But meanwhile you already made an error by asking help and an output.\n");
	}
}

void *Initiate_Session (unsigned int verbose)
{	/* Initialize the GMT Session and store the API pointer in a persistent variable */
	void *API = NULL;
	/* Initializing new GMT session with zero pad and a Matlab-acceptable replacement for the printf function */
	if ((API = GMT_Create_Session (MEX_PROG, 0U, (verbose << 2) + GMT_SESSION_NOEXIT + GMT_SESSION_EXTERNAL, GMTMEX_print_func)) == NULL)
		mexErrMsgTxt ("GMT: Failure to create new GMT session\n");

	if (!pPersistent) pPersistent = mxMalloc(sizeof(uintptr_t));
	pPersistent[0] = (uintptr_t)(API);
	mexMakeMemoryPersistent (pPersistent);

	return (API);
}

/* This is the function that is called when we type gmt in Matlab/Octave */
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int status = 0;                 /* Status code from GMT API */
	unsigned int first = 0;         /* Array ID of first command argument (not 0 when API-ID is first) */
	unsigned int verbose = 0;       /* Default verbose setting */
	unsigned int n_items = 0, pos = 0; /* Number of Matlab arguments (left and right) */
	size_t str_length = 0, k = 0;   /* Misc. counters */
	void *API = NULL;		/* GMT API control structure */
	struct GMT_OPTION *options = NULL; /* Linked list of module options */
	struct GMT_RESOURCE *X = NULL;  /* Array of information about Matlab args */
	char *cmd = NULL;               /* Pointer used to get the user's Matlab command */
	char *gtxt = NULL;              /* For debug printing of revised command */
	char *opt_args = NULL;		/* Pointer to the user's module options */
	const char *keys = NULL;	/* This module's option keys */
	char module[MODULE_LEN] = {""}; /* Name of GMT module to call */
	char name[GMT_STR16];           /* Name of GMT module to call */
	void *ptr = NULL;
	uintptr_t *pti = NULL;          /* To locally store the API address */

	/* 0. No arguments at all results in the GMT banner message */
	if (nrhs == 0) {
		usage (nlhs, nrhs);
		return;
	}

	/* 1. Check for the special commands create and help */
	
	if (nrhs == 1) {	/* This may be create or help */
		cmd = mxArrayToString (prhs[0]);
		if (!strncmp (cmd, "help", 4U) || !strncmp (cmd, "--help", 6U)) {
			usage (nlhs, 1);
			return;
		}
		if (!strncmp (cmd, "create", 6U)) {	/* Asked to create a new GMT session */
			if (nlhs > 1)	/* Asked for too much output, only 1 or 0 is allowed */
				mexErrMsgTxt ("GMT: Usage: gmt ('create') or API = gmt ('create');\n");
			if (pPersistent)                        /* See if have a GMT API pointer */
				API = (void *)pPersistent[0];
			if (API != NULL) {                      /* If another session still exists */
				mexPrintf ("GMT: A previous GMT session is still active. Ignoring your 'create' request.\n");
				if (nlhs) /* Return nothing */
					plhs[0] = mxCreateNumericMatrix (1, 0, mxUINT64_CLASS, mxREAL);
				return;
			}
			if ((gtxt = strstr (cmd, "-V"))) verbose = GMT_get_V (gtxt[2]);
			API = Initiate_Session (verbose);	/* Initializing a new GMT session */

			if (nlhs) {	/* Return the API adress as an integer (nlhs == 1 here) )*/
				plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
				pti = mxGetData(plhs[0]);
				*pti = *pPersistent;
			}

			mexAtExit(force_Destroy_Session);	/* Register an exit function. */
			return;
		}

		/* OK, neither create nor help, must be a single command with no arguments nor the API. So get it: */
		if (!pPersistent)
			mexErrMsgTxt ("GMT: You shouldn't have cleared this mex. Now the GMT5 session is lost (mem leaked).\n"); 
		API = (void *)pPersistent[0];	/* Get the GMT API pointer */
		if (API == NULL) mexErrMsgTxt ("GMT: This GMT5 session has already been destroyed, or currupted.\n"); 
	}
	else if (mxIsScalar_(prhs[0]) && mxIsUint64(prhs[0])) {
		/* Here, nrhs > 1 . If first arg is a scalar int, we assume it is the API memory address */
		pti = (uintptr_t *)mxGetData(prhs[0]);
		API = (void *)pti[0];	/* Get the GMT API pointer */
		first = 1;		/* Commandline args start at prhs[1] since prhs[0] had the API id argument */
	}
	else {		/* We still don't have the API, so we must get it from the past or initiate a new session */
		if (!pPersistent || (API = (void *)pPersistent[0]) == NULL)
			API = Initiate_Session (verbose);	/* Initializing new GMT session */
	}

	if (!cmd) 	/* First argument is the command string, e.g., 'blockmean -R0/5/0/5 -I1' or just 'destroy' */
		cmd = mxArrayToString (prhs[first]);

	if (!strncmp (cmd, "destroy", 7U)) {	/* Destroy the session */
		if (nlhs != 0)
			mexErrMsgTxt ("GMT: Usage is gmt ('destroy');\n");

		if (GMT_Destroy_Options (API, &options)) mexErrMsgTxt ("GMT: Failure to destroy GMT5 options\n");
		if (GMT_Destroy_Session (API)) mexErrMsgTxt ("GMT: Failure to destroy GMT5 session\n");
		*pPersistent = 0;	/* Wipe the persistent memory */
		return;
	}

	/* Here we have a GMT module call */
	
	/* 2. Get mex arguments, if any, and extract the GMT module name */
	str_length = strlen (cmd);				/* Length of command argument */
	for (k = 0; k < str_length && cmd[k] != ' '; k++);	/* Determine first space in command */
	if (k >= MODULE_LEN)
		mexErrMsgTxt ("GMT: Module name in command is too long\n");
	strncpy (module, cmd, k);				/* Isolate the module name in this string */

	/* 3. Convert mex command line arguments to a linked GMT option list */
	while (cmd[k] == ' ') k++;	/* Skip any spaces between module name and start of options */
	opt_args = (cmd[k]) ? &cmd[k] : NULL;
	if (opt_args && (options = GMT_Create_Options (API, 0, opt_args)) == NULL)
		mexErrMsgTxt ("GMT: Failure to parse GMT5 command options\n");

	/* 4. Preprocess to update GMT option lists and return info array X */
	if ((X = GMT_Encode_Options (API, module, ARG_MARKER, nlhs, &options, &n_items)) == NULL)
		mexErrMsgTxt ("GMT: Failure to encode mex command options\n");
	
	/* 5. Assign input (from mex) and output (from GMT) resources */
	
	for (k = 0; k < n_items; k++) {	/* Number of GMT containers involved in this module call */
		ptr = (X[k].direction == GMT_IN) ? (void *)prhs[X[k].pos+first+1] : (void *)plhs[X[k].pos];
		X[k].object = GMTMEX_Register_IO (API, X[k].family, X[k].geometry, X[k].direction, ptr, &X[k].object_ID);
		if (X[k].object == NULL || X[k].object_ID == GMT_NOTSET)
			mexErrMsgTxt("GMT: Failure to register the resource\n");
		if (GMT_Encode_ID (API, name, X[k].object_ID) != GMT_NOERROR) 	/* Make filename with embedded object ID */
			mexErrMsgTxt ("GMT: Failure to encode string\n");
		if (GMT_Expand_Option (API, X[k].option, ARG_MARKER, name) != GMT_NOERROR)	/* Replace ARG_MARKER in argument with name */
			mexErrMsgTxt ("GMT: Failure to expand filename marker\n");
	}
	
	/* 6. Run GMT module; give usage message if errors arise during parsing */
	gtxt = GMT_Create_Cmd (API, options);
	GMT_Report (API, GMT_MSG_DEBUG, "GMT_Encode_Options: Revised command after memory-substitution: %s\n", gtxt);
	GMT_Destroy_Cmd (API, &gtxt);	/* Only needed it for the above verbose */
	
	if ((status = GMT_Call_Module (API, module, GMT_MODULE_OPT, options)) != GMT_NOERROR)
		mexErrMsgTxt ("GMT: Module return with failure\n");

	/* 7. Hook up module GMT outputs to Matlab plhs array */
	
	for (k = 0; k < n_items; k++) {	/* Number of GMT containers involved in this module call */
		if (X[k].direction == GMT_OUT) {	/* Get results from GMT into Matlab arrays */
			if ((X[k].object = GMT_Retrieve_Data (API, X[k].object_ID)) == NULL)
				mexErrMsgTxt ("GMT: Error retrieving object from GMT\n");
			pos = X[k].pos;	/* Short-hand for index into the plhs[] array being returned to Matlab */
			switch (X[k].family) {	/* Determine what container we got */
				case GMT_IS_GRID:	/* A GMT grid; make it the pos'th output item */
					plhs[pos] = GMTMEX_Get_Grid (API, X[k].object);
					break;
				case GMT_IS_DATASET:	/* A GMT table; make it a matrix and the pos'th output item */
					plhs[pos] = GMTMEX_Get_Dataset (API, X[k].object);
					break;
				case GMT_IS_TEXTSET:	/* A GMT textset; make it a cell and the pos'th output item */
					plhs[pos] = GMTMEX_Get_Textset (API, X[k].object);
					break;
				case GMT_IS_CPT:	/* A GMT CPT; make it a colormap and the pos'th output item  */
					plhs[pos] = GMTMEX_Get_CPT (API, X[k].object);
					break;
				case GMT_IS_IMAGE:	/* A GMT Image; make it the pos'th output item  */
					plhs[pos] = GMTMEX_Get_Image (API, X[k].object);
					break;
				default:
					mexErrMsgTxt ("GMT: Unsupported data type\n");
					break;
			}
		}	/* else means we have an object used to pass arrays from Matlab into GMT */
		else {	/* Free any memory allocated outside of GMT */
			if (X[k].family == GMT_IS_TEXTSET)	/* Because Windows cannot stomach another DLL freeing strings */
				GMTMEX_Free_Textset (API, X[k].object);
		}
			
		if (GMT_Destroy_Data (API, &X[k].object) != GMT_NOERROR)
			mexErrMsgTxt ("GMT: Failed to destroy object used in the interface bewteen GMT and Matlab\n");
	}
	
	/* 8. Destroy linked option list */
	
	if (GMT_Destroy_Options (API, &options)) mexErrMsgTxt ("GMT: Failure to destroy GMT5 options\n");
	return;
}
