#	$Id$
#
#	Authors: Paul Wessel and Joaquim Luis
#	Updated: July 26, 2015

The Matlab/Octave API is built on top of the GMT C/C++ API
and is therefore decoupled from the GMT5 distribution.  We
anticipate this is how other APIs (Python, Julia) will be
done as well.  To build gmtmex you first must install
GMT >=5.2 from subversion or via an official release.

The single gmt.c program is linked with gmtmex_parser.c and
the GMT5 API library into a mex executable.  Parsing is
assisted via program-specific option information obtained via
the API function GMT_Encode_Options. A configure option will
select between Matlab versus Octave builds.

To build the gmt mex executable, take the following steps:

1. Run autoconf
2. Run configure --help and acquaint yourself with the various options.
3. Then, run configure with the selected options, in addition to using
   either --enable-matlab or --enable-octave (and --enable-octmode=mex|oct).
4. Finally, run make all and make install

To debug the mex please see the wiki for gmt-mex at
http://gmt.soest.hawaii.edu/projects/gmt-MATLAB-octave-api/wiki

GMT Setup:

Due to the problematic way MATLAB maintains its own old shared libraries
we must for now mess quite a bit to get things to work.  This is OK
for now:

1. Build and installed the OS X bundle for 5.2.
2. Run the set_osx_libs.sh; see instructions in the script.
3. Make a netcdf-3 version of gshhg using nccopy
4. Make sure config.mk lists the -lXgmt and correct path (2 above)
5. Run make to build gmt.mexmaci64
6. setenv DYLD_LIBRARY_PATH to the /Application/.../mex dir
7. Make changes to your gmt.conf file:
   a) Change default IO_NC4_CHUNK_SIZE to classic.
   b) Change DIR_GSHHG to the nc3 version
   c) Set GMT_CUSTOM_LIBS to /Application/.../mex/gmt/plugins/supplements.so
8. Start MATLAB from command line.

With those changes and limitations (cannot read netcdf4 grid) all
operations in GMT works except psconvert which calls gs and it
has newer shared libraries than MATLAB is willing to acknowledge.
Also, using DCW via pscoast does not work.
