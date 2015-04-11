#!/bin/bash
# $Id$
# Rebaptizing of shared libraries for OSX.  We will duplicate shared libraries
# needed by GMT and its dependencies that are distributed with Matlab.  We place
# these copies in a mex sub-directory where all GMT libs live.
#
# Set up for macports for now
# Must be run as sudo
SLIB=/usr/lib			# System library folder
GLIB=/usr/local/lib		# GMT library folder
PLIB=/opt/local/lib		# MacPort library folder
MLIB=${GLIB}/mex		# Folder with rebaptized libraries
#-----------------------------------------------------------------
mkdir -p ${MLIB}/gmt/plugins
HERE=`pwd`
cd $MLIB

# Z
cp ${PLIB}/libz.1.2.8.dylib libz-mex.1.2.8.dylib
ln -s libz-mex.1.2.8.dylib libz-mex.1.dylib
ln -s libz-mex.1.2.8.dylib libz-mex.dylib
chmod -h og+r libz-mex.1.dylib
chmod -h og+r libz-mex.dylib
install_name_tool -id libz-mex.dylib libz-mex.1.2.8.dylib
# CRYPTO
cp ${PLIB}/libcrypto.1.0.0.dylib libcrypto-mex.1.0.0.dylib
ln -s libcrypto-mex.1.0.0.dylib libcrypto-mex.dylib
chmod -h og+r libcrypto-mex.dylib
install_name_tool -id libcrypto-mex.dylib libcrypto-mex.1.0.0.dylib
install_name_tool -change ${PLIB}/libz.1.dylib ${MLIB}/libz-mex.1.dylib libcrypto-mex.1.0.0.dylib
# TIFF
cp ${PLIB}/libtiff.5.dylib libtiff-mex.5.dylib 
ln -s libtiff-mex.5.dylib libtiff-mex.dylib 
chmod -h og+r libtiff-mex.dylib 
install_name_tool -id libtiff-mex.dylib libtiff-mex.5.dylib
install_name_tool -change ${PLIB}/libz.1.dylib ${MLIB}/libz-mex.1.dylib libtiff-mex.5.dylib
# PROJ
cp ${PLIB}/libproj.9.dylib libproj-mex.9.dylib 
ln -s libproj-mex.9.dylib libproj-mex.dylib
chmod -h og+r libproj-mex.dylib
install_name_tool -id libproj-mex.dylib libproj-mex.9.dylib
# SSL
cp ${PLIB}/libssl.1.0.0.dylib libssl-mex.1.0.0.dylib 
ln -s libssl-mex.1.0.0.dylib libssl-mex.dylib
chmod -h og+r libssl-mex.dylib
install_name_tool -id libssl-mex.dylib libssl-mex.1.0.0.dylib
install_name_tool -change ${PLIB}/libcrypto.1.0.0.dylib ${MLIB}/libcrypto-mex.1.0.0.dylib libssl-mex.1.0.0.dylib
install_name_tool -change ${PLIB}/libz.1.dylib		${MLIB}/libz-mex.1.dylib	  libssl-mex.1.0.0.dylib
# GEOTIFF
cp ${PLIB}/libgeotiff.2.dylib libgeotiff-mex.2.dylib 
ln -s libgeotiff-mex.2.dylib libgeotiff-mex.dylib 
chmod -h og+r libgeotiff-mex.dylib 
install_name_tool -id libgeotiff-mex.dylib libgeotiff-mex.2.dylib
install_name_tool -change ${PLIB}/libproj.9.dylib       ${MLIB}/libproj-mex.9.dylib       libgeotiff-mex.2.dylib
install_name_tool -change ${PLIB}/libtiff.5.dylib       ${MLIB}/libtiff-mex.5.dylib       libgeotiff-mex.2.dylib
install_name_tool -change ${PLIB}/libz.1.dylib 		${MLIB}/libz-mex.1.dylib	  libgeotiff-mex.2.dylib
# CURL
cp ${PLIB}/libcurl.4.dylib libcurl-mex.4.dylib 
ln -s libcurl-mex.4.dylib libcurl-mex.dylib
chmod -h og+r libcurl-mex.dylib
install_name_tool -id libcurl-mex.dylib libcurl-mex.4.dylib
install_name_tool -change ${PLIB}/libcrypto.1.0.0.dylib ${MLIB}/libcrypto-mex.1.0.0.dylib libcurl-mex.4.dylib
install_name_tool -change ${PLIB}/libz.1.dylib ${MLIB}/libz-mex.1.dylib libcurl-mex.4.dylib
# HDF5
cp ${PLIB}/libhdf5.9.dylib libhdf5-mex.9.dylib 
ln -s libhdf5-mex.9.dylib libhdf5-mex.dylib
chmod -h og+r libhdf5-mex.dylib
install_name_tool -id libhdf5-mex.dylib libhdf5-mex.9.dylib
install_name_tool -change ${PLIB}/libz.1.dylib ${MLIB}/libz-mex.1.dylib libhdf5-mex.9.dylib
# HDF5_HL
cp ${PLIB}/libhdf5_hl.9.dylib libhdf5_hl-mex.9.dylib 
ln -s libhdf5_hl-mex.9.dylib libhdf5_hl-mex.dylib
chmod -h og+r libhdf5_hl-mex.dylib
install_name_tool -id libhdf5_hl-mex.dylib libhdf5_hl-mex.9.dylib
install_name_tool -change ${PLIB}/libhdf5.9.dylib       ${MLIB}/libhdf5-mex.9.dylib       libhdf5_hl-mex.9.dylib
install_name_tool -change ${PLIB}/libz.1.dylib		${MLIB}/libz-mex.1.dylib	  libhdf5_hl-mex.9.dylib
# NETCDF:
cp ${PLIB}/libnetcdf.7.dylib libnetcdf-mex.7.dylib 
ln -s libnetcdf-mex.7.dylib libnetcdf-mex.dylib
chmod -h og+r libnetcdf-mex.dylib
install_name_tool -id libnetcdf-mex.dylib libnetcdf-mex.7.dylib
install_name_tool -change ${PLIB}/libhdf5.9.dylib       ${MLIB}/libhdf5-mex.9.dylib       libnetcdf-mex.7.dylib
install_name_tool -change ${PLIB}/libhdf5_hl.9.dylib    ${MLIB}/libhdf5_hl-mex.9.dylib    libnetcdf-mex.7.dylib
install_name_tool -change ${PLIB}/libcurl.4.dylib       ${MLIB}/libcurl-mex.4.dylib       libnetcdf-mex.7.dylib
install_name_tool -change ${PLIB}/libz.1.dylib		${MLIB}/libz-mex.1.dylib	  libnetcdf-mex.7.dylib
# GDAL:
cp ${PLIB}/libgdal.1.dylib libgdal-mex.1.dylib
ln -s libgdal-mex.1.dylib libgdal-mex.dylib
chmod -h og+r libgdal-mex.dylib
install_name_tool -id libgdal-mex.dylib libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libnetcdf.7.dylib     ${MLIB}/libnetcdf-mex.7.dylib     libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libhdf5.9.dylib       ${MLIB}/libhdf5-mex.9.dylib       libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libgeotiff.2.dylib    ${MLIB}/libgeotiff-mex.2.dylib    libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libcurl.4.dylib       ${MLIB}/libcurl-mex.4.dylib       libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libssl.1.0.0.dylib    ${MLIB}/libssl-mex.1.0.0.dylib    libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libproj.9.dylib       ${MLIB}/libproj-mex.9.dylib       libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libgeotiff.2.dylib    ${MLIB}/libgeotiff-mex.2.dylib    libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libtiff.5.dylib       ${MLIB}/libtiff-mex.5.dylib       libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libcrypto.1.0.0.dylib ${MLIB}/libcrypto-mex.1.0.0.dylib libgdal-mex.1.dylib
install_name_tool -change ${PLIB}/libz.1.dylib		${MLIB}/libz-mex.1.dylib	  libgdal-mex.1.dylib
# PSL:
cp ${GLIB}/libpsl.5.2.0.dylib libpsl-mex.5.2.0.dylib
ln -s libpsl-mex.5.2.0.dylib libpsl-mex.5.dylib
chmod -h og+r libpsl-mex.5.dylib
ln -s libpsl-mex.5.dylib     libpsl-mex.dylib
chmod -h og+r libpsl-mex.dylib
install_name_tool -id libpsl-mex.dylib libpsl-mex.5.2.0.dylib
install_name_tool -change ${SLIB}/libz.1.dylib	    ${MLIB}/libz-mex.1.dylib	  libpsl-mex.5.2.0.dylib
# GMT:
cp ${GLIB}/libgmt.5.2.0.dylib libgmt-mex.5.2.0.dylib
ln -s libgmt-mex.5.2.0.dylib libgmt-mex.5.dylib
chmod -h og+r libgmt-mex.5.dylib
ln -s libgmt-mex.5.dylib     libgmt-mex.dylib
chmod -h og+r libgmt-mex.dylib
install_name_tool -id libgmt-mex.dylib libgmt-mex.5.2.0.dylib
install_name_tool -change ${PLIB}/libnetcdf.7.dylib ${MLIB}/libnetcdf-mex.7.dylib libgmt-mex.5.2.0.dylib
install_name_tool -change ${PLIB}/libgdal.1.dylib   ${MLIB}/libgdal-mex.1.dylib   libgmt-mex.5.2.0.dylib
install_name_tool -change ${PLIB}/libgdal.1.dylib   ${MLIB}/libgdal-mex.1.dylib   libgmt-mex.5.2.0.dylib
install_name_tool -change ${SLIB}/libz.1.dylib	    ${MLIB}/libz-mex.1.dylib	  libgmt-mex.5.2.0.dylib
install_name_tool -change ${GLIB}/libpsl.5.dylib    ${MLIB}/libpsl-mex.5.dylib	  libgmt-mex.5.2.0.dylib
cp ${GLIB}/gmt/plugins/supplements.so ${MLIB}/gmt/plugins
cd ${MLIB}/gmt/plugins
install_name_tool -change ${GLIB}/libgmt.5.dylib    ${MLIB}/libgmt-mex.5.dylib    supplements.so
install_name_tool -change ${GLIB}/libpsl.5.dylib    ${MLIB}/libpsl-mex.5.dylib    supplements.so
install_name_tool -change ${PLIB}/libnetcdf.7.dylib ${MLIB}/libnetcdf-mex.7.dylib supplements.so
install_name_tool -change ${PLIB}/libgdal.1.dylib   ${MLIB}/libgdal-mex.1.dylib   supplements.so
install_name_tool -change ${SLIB}/libz.1.dylib	    ${MLIB}/libz-mex.1.dylib	  supplements.so
cd ${MLIB}
cd ..
chmod -R og+r mex
cd $HERE
echo "Done"
