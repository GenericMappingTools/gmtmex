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
   either --enable-matlab or --enable-octave (and --enable-octmode=mex|oct).
4. Finally, run make all and make install

To debug the mex please see the wiki for gmt-mex at
http://gmt.soest.hawaii.edu/projects/gmt-matlab-octave-api/wiki

GMT Setup:

Due to the problematic way matlab maintains its own old shared libraries
we must for now mess quite a bit to get things to work.  This is OK
for now:

1. Build and installed the OS X bundle for 5.2.
2. Run the prep_osxbundle.sh; see instructions in the script.
3. Make a netcdf-3 version of gshhg usinc nccopy
4. Change default IO_NC4_CHUNK_SIZE to classic.
5. Change DIR_GSHHG to the nc3 version
6. Make sure config.mk lists the -lXgmt and correct path (2 above)
7. Run make to build gmt.mexmaci64
8. setenv DYLD_LIBRARY_PATH to the /Application/.../mex dir
9. Start matlab from command line.

With those changes and limitations (cannot read netcdf4 grid) all
operations in GMT works except psconvert which calls gs and it
has newer shared libraries than matlab is willing to acknowledge.

