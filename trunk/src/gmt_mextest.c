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
	int nlhs = 0, nrhs = 0, k;				/* Simulated counts */
	int start, quotes;
	size_t str_length;				/* Misc. counters */
	struct GMTAPI_CTRL *API = NULL;			/* GMT API control structure */
	struct GMT_OPTION *options = NULL;		/* Linked list of options */
	struct GMTMEX *X = NULL;			/* Array of information about Matlab args */
	char *str = NULL;				/* Pointer used to get Matlab command */
	char module[BUFSIZ] = {""}, cmd[BUFSIZ] = {""};	/* Name of GMT module to call and the command */
	mxArray *plhs[5] = {NULL, NULL, NULL, NULL, NULL};	/* Simulated pointers to Matlab arrays */
	const mxArray *prhs[5] = {NULL, NULL, NULL, NULL, NULL};

	if (argc != 2) {
		fprintf (stderr, "\ngmt_mextest - Test mex argument parsing\n\n");
		fprintf (stderr, "usage: Single double-quoted GMT/MEX command of the forms\n");
		fprintf (stderr, "	\"gmt ('module -options');\"\n");
		fprintf (stderr, "	\"gmt ('module -options', X);\"\n");
		fprintf (stderr, "	\"A = gmt ('module -options');\"\n");
		fprintf (stderr, "	\"A = gmt ('module -options', X);\"\n");
		fprintf (stderr, "	\"[A, B, ...] = gmt ('module -options', X, Y, ...);\"\n");
		exit (-1);
	}

	/* 1. Initializing new GMT session with zero pad */
	if ((API = GMT_Create_Session ("GMT5", 0U, 3U, NULL)) == NULL) fprintf (stderr, "Failure to create GMT5 Session\n");

	/* Expect a single argument (due to enclosing double quotes) of these forms:
		gmt ('module -options');
		gmt ('module -options', X)
	 	A = gmt ('module -options');
	 	A = gmt ('module -options', X);
		[A,...] = gmt ('module -options', X, ...);
	 */
	
	if (strchr (argv[1], '=')) {	/* Gave output arguments via <left> = gmt ('   '); */
		nlhs = 1;	/* At least one output argument given */
		if (argv[1][0] == '[') {	/* Gave several output arguments in [] */
			k = 1;
			while (argv[1][k] && argv[1][k] != ']') {
				if (argv[1][k] == ',' || (argv[1][k] == ' ' && argv[1][k-1] != ',')) nlhs++;
				k++;
			}
		}
	}
	quotes = start = k = 0;
	while (quotes < 2 && argv[1][k]) {	/* Wind to 2nd single quote */
		if (start == 0 && argv[1][k] == '\'') start = k + 1;
		if (argv[1][k++] == '\'') quotes++;
	}
	strncpy (cmd, &argv[1][start], k-start-1);
	if (argv[1][k] != ')') {	/* Input arguments given */
		nrhs = 0;	/* At least one, count commas to determine total input items */
		while (argv[1][k] && argv[1][k] != ')') {
			if (argv[1][k++] == ',') nrhs++;
		}
	}

	/* 2. Get mex arguments, if any, and extract the GMT module name */
	str_length = strlen  (cmd);				/* Length of command argument */
	for (k = 0; k < str_length && cmd[k] != ' '; k++);	/* Determine first space in command */
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
	
	/* Print out expanded command line */
	if (nlhs) {
		if (nlhs == 1)
			printf ("A");
		else {
			printf ("[");
			for (k = 0; k < nlhs; k++) {
				if (k) putchar (' ');
				printf ("%c", 'A' + k);
			}
			printf ("]");
		}
		printf (" = ");
	}
	printf ("gmt (%s); %% Arguments: %d output %d input\n", revised_cmd, nlhs, nrhs);
	
	/* 6. Fake Run GMT module; g */

	/* 7. Hook up module output to Matlab plhs arguments */
	if (GMTMEX_post_process (API, X, n_items, plhs)) fprintf (stderr, "Failure to extract GMT5-produced data\n");
	
	/* 8. Destroy linked option list */
	if (GMT_Destroy_Options (API, &options)) fprintf (stderr, "Failure to destroy GMT5 options\n");

	/* 9. Destroy GMT session */
	if (GMT_Destroy_Session (API)) fprintf (stderr, "Failure to destroy GMT5 session\n");
	
	exit (0);
}
