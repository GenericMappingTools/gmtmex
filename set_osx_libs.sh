#!/bin/bash
#	$Id$
# Accept name of GMT distro in /Application, then
# 1) rebaptizes the libs
# 2) sets up links from /usr/local
# 3) prints the DYLD_LIBRARY_PATH you need to set
#
# E.g., (from gmt-mex/trunk): set_osx_libs.sh 
HERE=`pwd`
DIR="/Applications/$1"
cd $DIR/Contents/Resources/lib
$HERE/prep_osxbundle.sh
cd /usr/local/include
sudo rm -rf gmt
sudo ln -s $DIR/Contents/Resources/include/gmt gmt
sudo chmod -h og+r gmt
cd ../lib
sudo rm -f libgmt.dylib libpsl.dylib
sudo ln -s $DIR/Contents/Resources/lib/mex/libXgmt.dylib libgmt.dylib
sudo ln -s $DIR/Contents/Resources/lib/mex/libXpsl.dylib libpsl.dylib
sudo chmod -h og+r libgmt.dylib libpsl.dylib
echo "export DYLD_LIBRARY_PATH=$DIR/Contents/Resources/lib/mex"
