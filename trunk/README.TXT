#	$Id$
#
The new Matlab/Octave API will be built on top of the GMT C API only.  Thus, we
have decoupled its development from the GMT5 distribution.  We anticipate this
is how the other APIs (Python, Julia) will be done as well.  To build gmtmex
you first must update to the very latest GMT5 version from subversion and
build it.

The plan is to develop a generic GMTMEX_parser that is called by a dummy program
called mexprogram.c (modelled after gmtprogram.c).  Scripts will loop over all
GMT program names and supply their names to the process.  Program-specific option
information is kept in the mexproginfo.txt and we use this to build the API. Once
we get this right the entire API (i.e., mex files for each GMT program) will be
built automatically.  A configure option will select Matlab versus Octave builds.

For testing of parser, see gmt_mextest which you can run like this:

gmt_mextest nlhs nrhs "mex command line"

e.g.

gmt_mextest 1 1 "surface -R0/16/0/16 -I1 -T0.25 -V"

and it should determine that matlab provides input array and wants output grid returned.
The output should be

'surface -R0/16/0/16 -I1 -T0.25 -V @GMTAPI@-000000 -G@GMTAPI@-000001'
