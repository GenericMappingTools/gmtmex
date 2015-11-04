#	$Id$
#
#	Authors: Paul Wessel and Joaquim Luis
#	Updated: Nov 3, 2015

The Matlab/Octave API is built on top of the GMT C/C++ API
and is therefore decoupled from the GMT5 distribution.  We
anticipate this is how other APIs (Python, Julia) will be
done as well.  To build gmtmex you first must install
GMT >=5.2 from subversion or via an official release.

The single gmtmex.c program is linked with gmtmex_parser.c and
the GMT5 API library into a mex executable. A configure option
selects between Matlab versus Octave builds.  The script
gmt.m is responsible for calling the gmtmex executable.

GMT Setup:

Due to the problematic way MATLAB maintains its own old shared
libraries (such as netcdf, hdf5, and others) we must for now
apply a messy workaround to get things to work.  Basically, we
duplicate the bundle installation to /opt/gmt and change the
names of all the shared libraries.  The steps below are thus
required:

1. Build (or more likely just install) the OS X bundle for 5.2.
2. Run the gmt_prepmex.sh script in the share/tools folder.
   This duplicates the gmt installation to /opt/gmt.
3. Run gmtswitch and select /opt/gmt as the active version.
4. Run autoconf and configure --enable-matlab.
5. Run make to build gmtmex.mexmaci64
6. Make changes to your gmt.conf file:
     Set GMT_CUSTOM_LIBS to /Applications/GMT-5.2.0xxxxx.app/Contents/Resources/lib/gmt/plugins/
7. Adjust your MATLAB's path so that the gmt.m and gmtmex are found.

To debug the mex please see the wiki for gmt-mex at
http://gmt.soest.hawaii.edu/projects/gmt-MATLAB-octave-api/wiki

