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
 * This is a test program for the testing the completion of the implicit options. 
 *
 * Version:	5.2
 * Created:	9-Feb-2015
 *
 */

#include "gmt.h"
#include <string.h>

#define ARG_MARKER	'$'	/* Character that indicates an memory reference to data */

static const char *GMT_family[] = {"Data Table", "Text Table", "GMT Grid", "CPT Table", "GMT Image", "GMT Vector", "GMT Matrix", "GMT Coord"};
static const char *GMT_geometry[] = {"Not Set", "Point", "Line", "Polygon", "Point|Line|Poly", "Line|Poly", "Surface", "Non-Geographical"};
static const char *GMT_direction[] = {"Input", "Output"};
unsigned int gmtry (unsigned int geometry) {
	/* Return index to text representation in GMT_geometry[] */
	if (geometry == GMT_IS_POINT)   return 1;
	if (geometry == GMT_IS_LINE)    return 2;
	if (geometry == GMT_IS_POLY)    return 3;
	if (geometry == GMT_IS_PLP)     return 4;
	if ((geometry & GMT_IS_LINE) && (geometry & GMT_IS_POLY)) return 5;
	if (geometry == GMT_IS_SURFACE) return 6;
	if (geometry == GMT_IS_NONE)    return 7;
	return 0;
}

int main (int argc, char *argv[]) {
	unsigned int g;
	unsigned int n_items = 0;	/* Number of Matlab arguments (left and right) */
	unsigned int k;			/* Simulated counts */
	int start, quotes;
	size_t str_length;				/* Misc. counters */
	void *API = NULL;				/* GMT API control structure */
	struct GMT_OPTION *options = NULL;		/* Linked list of options */
	struct GMT_RESOURCE *X = NULL;			/* Array of information about Matlab args */
	char *str = NULL;				/* Pointer used to get Matlab command */
	char module[BUFSIZ] = {""}, cmd[BUFSIZ] = {""};	/* Name of GMT module to call and the command */
	char revised_cmd[BUFSIZ];			/* Show revised command when testing only */

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
	if ((API = GMT_Create_Session ("GMT5", 0, (GMT_MSG_DEBUG << 2)+GMT_SESSION_NOEXIT+GMT_SESSION_EXTERNAL, NULL)) == NULL) fprintf (stderr, "Failure to create GMT5 Session\n");

	/* Expect a single argument (due to enclosing double quotes) of these forms:
		gmt ('module -options');
		gmt ('module -options', X)
	 	A = gmt ('module -options');
	 	A = gmt ('module -options', X);
		[A,...] = gmt ('module -options', X, ...);
	 */
	
	quotes = start = k = 0;
	while (quotes < 2 && argv[1][k]) {	/* Wind to 2nd single quote */
		if (start == 0 && argv[1][k] == '\'') start = k + 1;
		if (argv[1][k++] == '\'') quotes++;
	}
	strncpy (cmd, &argv[1][start], k-start-1);

	/* 2. Get mex arguments, if any, and extract the GMT module name */
	str_length = strlen  (cmd);				/* Length of command argument */
	for (k = 0; k < str_length && cmd[k] != ' '; k++);	/* Determine first space in command */
	strncpy (module, cmd, k);				/* Isolate the module name in this string */

	/* 3. Convert mex command line arguments to a linked option list */
	str = (k < str_length) ? &cmd[k+1] : NULL;
	if ((options = GMT_Create_Options (API, 0, str)) == NULL)
		fprintf (stderr, "Failure to parse GMT5 command options\n");

	/* 4. Preprocess and update GMT option lists, and return X info array */
	if ((X = GMT_Encode_Options (API, module, ARG_MARKER, &options, &n_items)) == NULL)
		fprintf (stderr, "Failure to encode mex command options\n");

	printf ("Revised command: %s\n", revised_cmd);
	
	for (k = 0; k < n_items; k++) {
		g = gmtry (X[k].geometry);
		fprintf (stderr, "GMT_RESOURCE item %d:", k);
		fprintf (stderr, " Id = %d", X[k].object_ID);
		fprintf (stderr, " %s", GMT_direction[X[k].direction]);
		fprintf (stderr, " %s", GMT_family[X[k].family]);
		fprintf (stderr, " [%s]", GMT_geometry[g]);
		fprintf (stderr, " in position: %d\n", X[k].pos);
	}

	free (X);
	/* 8. Destroy linked option list */
	if (GMT_Destroy_Options (API, &options)) fprintf (stderr, "Failure to destroy GMT5 options\n");

	/* 9. Destroy GMT session */
	if (GMT_Destroy_Session (API)) fprintf (stderr, "Failure to destroy GMT5 session\n");
	
	exit (0);
}
