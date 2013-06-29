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
 * This is a test program for the GMTMEX_pre_process function. No
 * Matlab/Octave libraries or actions are involved 
 *
 * Version:	5
 * Created:	28-Jun-2013
 *
 */

#include "gmtmex.h"

int main (int argc, char *argv[]) {
	int module_id;					/* Module ID */
	int n_items = 0;				/* Number of Matlab arguments (left and right) */
	int nlhs, nrhs;					/* Simulated counts */
	size_t str_length, k;				/* Misc. counters */
	struct GMTAPI_CTRL *API = NULL;			/* GMT API control structure */
	struct GMT_OPTION *options = NULL;		/* Linked list of options */
	struct GMTMEX *X = NULL;			/* Array of information about Matlab args */
	char *cmd = NULL, *str = NULL;				/* Pointer used to get Matlab command */
	char module[BUFSIZ];				/* Name of GMT module to call */
	mxArray *plhs[5] = {NULL, NULL, NULL, NULL, NULL};	/* Simulated pointers to Matlab arrays */
	const mxArray *prhs[5] = {NULL, NULL, NULL, NULL, NULL};

	if (argc != 4) {
		fprintf (stderr, "\ngmt_mextest - Test mex argument parsing\n");
		fprintf (stderr, "usage: nlhs nrhs \"mex-command\"\n");
		exit (-1);
	}

	/* 1. Initializing new GMT session with zero pad */
	if ((API = GMT_Create_Session ("GMT5", 0U, 1U, NULL)) == NULL) fprintf (stderr, "Failure to create GMT5 Session\n");

	nlhs = atoi (argv[1]);
	nrhs = atoi (argv[2]);
	cmd = argv[3];
	
	/* 2. Get mex arguments, if any, and extract the GMT module name */
	str_length = strlen  (cmd);				/* Length of command argument */
	for (k = 0; k < str_length && cmd[k] != ' '; k++);	/* Determine first space in command */
	memset ((void *)module, 0, BUFSIZ*sizeof (char));	/* Initialize module name to blank */
	strncpy (module, cmd, k);				/* Isolate the module name in this string */

	/* 3. Determine the GMT module ID, or list module usages and return if module is not found */
	if ((module_id = GMTMEX_find_module (API, module)) == -1) {
		GMT_Call_Module (API, NULL, GMT_MODULE_PURPOSE, NULL);
		if (GMT_Destroy_Session (API)) fprintf (stderr, "Failure to destroy GMT5 session\n");
		exit (-1);
	}

	/* 4. Convert mex command line arguments to a linked option list */
	str = (k < str_length) ? &cmd[k+1] : NULL;
	if ((options = GMT_Create_Options (API, 0, str)) == NULL)
		fprintf (stderr, "Failure to parse GMT5 command options\n");

	/* 5. Parse the mex command, update GMT option lists, and register in/out resources, and return X array */
	if ((n_items = GMTMEX_pre_process (API, module, plhs, nlhs, prhs, nrhs-1, keys[module_id], &options, &X)) < 0)
		fprintf (stderr, "Failure to parse mex command options\n");
	
	/* 6. Fake Run GMT module; g */

	/* 7. Hook up module output to Matlab plhs arguments */
	if (GMTMEX_post_process (API, X, n_items, plhs)) fprintf (stderr, "Failure to extract GMT5-produced data\n");
	
	/* 8. Destroy linked option list */
	if (GMT_Destroy_Options (API, &options)) fprintf (stderr, "Failure to destroy GMT5 options\n");

	/* 9. Destroy GMT session */
	if (GMT_Destroy_Session (API)) fprintf (stderr, "Failure to destroy GMT5 session\n");
	
	exit (0);
}
