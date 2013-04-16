#include "gmt.h"
int GMTMEX_parser (void *API, void *plhs[], int nlhs, void *prhs[], int nrhs, char *keys, struct GMT_OPTION *head);

#include "gmtmex_id.h"
#include "gmtmex_keys.h"
#include "gmtmex_progs.h"

int main (int argc, char *argv[])
{	/* Test driver for GMTMEX_parser: Use it like this:
	 * mexparser k_program 'typical mex args' */

	enum GMT_prog_enum id = k_dummy;
	void *API = NULL;			/* GMT API control structure */
	struct GMT_OPTION *options = NULL;	/* Linked list of options */
	char *cmd = NULL;

	id = atoi (argv[1]);	/* ID number of program to test */

	/* 1. Initializing new GMT session with no grid pad as default */
	if ((API = GMT_Create_Session ("GMT/MEX-API", 0U, 0U)) == NULL) fprintf (stderr, "Failure to create GMT Session\n");

	/* 2. Convert command line arguments to local linked option list */
	if ((options = GMT_Create_Options (API, 0, argv[2])) == NULL) fprintf (stderr, "Failure to parse GMT command options\n");

	/* 3. Parse the mex command, update GMT option lists, register in/out resources */
	if (GMTMEX_parser (API, NULL, 0, NULL, 0, keys[id], options)) fprintf (stderr, "Failure to parse mex command options\n");
	
	cmd = GMT_Create_Cmd (API, options);
	
	fprintf (stderr, "Call %s (API, -1, \"%s\")\n", func[id], cmd);
	free ((void *)cmd);
	
	/* 4. Destroy local linked option list */
	if (GMT_Destroy_Options (API, &options)) fprintf (stderr, "Failure to destroy GMT options\n");

	/* 5. Destroy GMT session */
	if (GMT_Destroy_Session (API)) fprintf (stderr, "Failure to destroy GMT session\n");

	exit (EXIT_SUCCESS);		/* Return the status from mex FUNC */
}
