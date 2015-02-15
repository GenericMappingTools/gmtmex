#	$Id$
#
#	Authors: Paul Wessel and Joaquim Luis
#	Updated: Feb 14, 2015

The Matlab/Octave API is built on top of the GMT C/C++ API
and is therefore decoupled from the GMT5 distribution.  We
anticipate this is how other APIs (Python, Julia) will be
done as well.  To build gmtmex you first must install
GMT >=5.2 from subversion or via an official release.

The single gmt.c program is linked with gmtmex_parser.c and
the GMT5 API library into a mex executable.  Parsing is
assisted via program-specific option information obtained via
GMT_Get_Moduleinfo. A configure option will select Matlab
versus Octave builds.

To build the gmt mex executable, take the following steps:

1. Run autoconf
2. Run configure --help and acquaint yourself with the various options.
3. Then, run configure with the selected options, in addition to using
   either --enable-matlab or --enable-debug.
4. Finally, run make all and make install

To debug the mex please see the wiki for gmt-mex at
http://gmt.soest.hawaii.edu/projects/gmt-matlab-octave-api/wiki
