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
cp ${PLIB}/libz.1.2.8.dylib libz_mex.1.2.8.dylib
ln -s libz_mex.1.2.8.dylib libz_mex.1.dylib
ln -s libz_mex.1.2.8.dylib libz_mex.dylib
chmod -h og+r libz_mex.1.dylib
chmod -h og+r libz_mex.dylib
install_name_tool -id libz_mex.dylib libz_mex.1.2.8.dylib
# CRYPTO
cp ${PLIB}/libcrypto.1.0.0.dylib libcrypto_mex.1.0.0.dylib
ln -s libcrypto_mex.1.0.0.dylib libcrypto_mex.dylib
chmod -h og+r libcrypto_mex.dylib
install_name_tool -id libcrypto_mex.dylib libcrypto_mex.1.0.0.dylib
install_name_tool -change ${PLIB}/libz.1.dylib ${MLIB}/libz_mex.1.dylib libcrypto_mex.1.0.0.dylib
# TIFF
cp ${PLIB}/libtiff.5.dylib libtiff_mex.5.dylib 
ln -s libtiff_mex.5.dylib libtiff_mex.dylib 
chmod -h og+r libtiff_mex.dylib 
install_name_tool -id libtiff_mex.dylib libtiff_mex.5.dylib
install_name_tool -change ${PLIB}/libz.1.dylib ${MLIB}/libz_mex.1.dylib libtiff_mex.5.dylib
# PROJ
cp ${PLIB}/libproj.9.dylib libproj_mex.9.dylib 
ln -s libproj_mex.9.dylib libproj_mex.dylib
chmod -h og+r libproj_mex.dylib
install_name_tool -id libproj_mex.dylib libproj_mex.9.dylib
# SSL
cp ${PLIB}/libssl.1.0.0.dylib libssl_mex.1.0.0.dylib 
ln -s libssl_mex.1.0.0.dylib libssl_mex.dylib
chmod -h og+r libssl_mex.dylib
install_name_tool -id libssl_mex.dylib libssl_mex.1.0.0.dylib
install_name_tool -change ${PLIB}/libcrypto.1.0.0.dylib ${MLIB}/libcrypto_mex.1.0.0.dylib libssl_mex.1.0.0.dylib
install_name_tool -change ${PLIB}/libz.1.dylib		${MLIB}/libz_mex.1.dylib	  libssl_mex.1.0.0.dylib
# GEOTIFF
cp ${PLIB}/libgeotiff.2.dylib libgeotiff_mex.2.dylib 
ln -s libgeotiff_mex.2.dylib libgeotiff_mex.dylib 
chmod -h og+r libgeotiff_mex.dylib 
install_name_tool -id libgeotiff_mex.dylib libgeotiff_mex.2.dylib
install_name_tool -change ${PLIB}/libproj.9.dylib       ${MLIB}/libproj_mex.9.dylib       libgeotiff_mex.2.dylib
install_name_tool -change ${PLIB}/libtiff.5.dylib       ${MLIB}/libtiff_mex.5.dylib       libgeotiff_mex.2.dylib
install_name_tool -change ${PLIB}/libz.1.dylib 		${MLIB}/libz_mex.1.dylib	  libgeotiff_mex.2.dylib
# CURL
cp ${PLIB}/libcurl.4.dylib libcurl_mex.4.dylib 
ln -s libcurl_mex.4.dylib libcurl_mex.dylib
chmod -h og+r libcurl_mex.dylib
install_name_tool -id libcurl_mex.dylib libcurl_mex.4.dylib
install_name_tool -change ${PLIB}/libcrypto.1.0.0.dylib ${MLIB}/libcrypto_mex.1.0.0.dylib libcurl_mex.4.dylib
install_name_tool -change ${PLIB}/libz.1.dylib ${MLIB}/libz_mex.1.dylib libcurl_mex.4.dylib
# HDF5
cp ${PLIB}/libhdf5.9.dylib libhdf5_mex.9.dylib 
ln -s libhdf5_mex.9.dylib libhdf5_mex.dylib
chmod -h og+r libhdf5_mex.dylib
install_name_tool -id libhdf5_mex.dylib libhdf5_mex.9.dylib
install_name_tool -change ${PLIB}/libz.1.dylib ${MLIB}/libz_mex.1.dylib libhdf5_mex.9.dylib
# HDF5_HL
cp ${PLIB}/libhdf5_hl.9.dylib libhdf5_hl_mex.9.dylib 
ln -s libhdf5_hl_mex.9.dylib libhdf5_hl_mex.dylib
chmod -h og+r libhdf5_hl_mex.dylib
install_name_tool -id libhdf5_hl_mex.dylib libhdf5_hl_mex.9.dylib
install_name_tool -change ${PLIB}/libhdf5.9.dylib       ${MLIB}/libhdf5_mex.9.dylib       libhdf5_hl_mex.9.dylib
install_name_tool -change ${PLIB}/libz.1.dylib		${MLIB}/libz_mex.1.dylib	  libhdf5_hl_mex.9.dylib
# NETCDF:
cp ${PLIB}/libnetcdf.7.dylib libnetcdf_mex.7.dylib 
ln -s libnetcdf_mex.7.dylib libnetcdf_mex.dylib
chmod -h og+r libnetcdf_mex.dylib
install_name_tool -id libnetcdf_mex.dylib libnetcdf_mex.7.dylib
install_name_tool -change ${PLIB}/libhdf5.9.dylib       ${MLIB}/libhdf5_mex.9.dylib       libnetcdf_mex.7.dylib
install_name_tool -change ${PLIB}/libhdf5_hl.9.dylib    ${MLIB}/libhdf5_hl_mex.9.dylib    libnetcdf_mex.7.dylib
install_name_tool -change ${PLIB}/libcurl.4.dylib       ${MLIB}/libcurl_mex.4.dylib       libnetcdf_mex.7.dylib
install_name_tool -change ${PLIB}/libz.1.dylib		${MLIB}/libz_mex.1.dylib	  libnetcdf_mex.7.dylib
# GDAL:
cp ${PLIB}/libgdal.1.dylib libgdal_mex.1.dylib
ln -s libgdal_mex.1.dylib libgdal_mex.dylib
chmod -h og+r libgdal_mex.dylib
install_name_tool -id libgdal_mex.dylib libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libnetcdf.7.dylib     ${MLIB}/libnetcdf_mex.7.dylib     libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libhdf5.9.dylib       ${MLIB}/libhdf5_mex.9.dylib       libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libgeotiff.2.dylib    ${MLIB}/libgeotiff_mex.2.dylib    libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libcurl.4.dylib       ${MLIB}/libcurl_mex.4.dylib       libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libssl.1.0.0.dylib    ${MLIB}/libssl_mex.1.0.0.dylib    libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libproj.9.dylib       ${MLIB}/libproj_mex.9.dylib       libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libgeotiff.2.dylib    ${MLIB}/libgeotiff_mex.2.dylib    libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libtiff.5.dylib       ${MLIB}/libtiff_mex.5.dylib       libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libcrypto.1.0.0.dylib ${MLIB}/libcrypto_mex.1.0.0.dylib libgdal_mex.1.dylib
install_name_tool -change ${PLIB}/libz.1.dylib		${MLIB}/libz_mex.1.dylib	  libgdal_mex.1.dylib
# PSL:
cp ${GLIB}/libpostscriptlight.5.2.0.dylib libpostscriptlight_mex.5.2.0.dylib
ln -s libpostscriptlight_mex.5.2.0.dylib libpostscriptlight_mex.5.dylib
chmod -h og+r libpostscriptlight_mex.5.dylib
ln -s libpostscriptlight_mex.5.dylib     libpostscriptlight_mex.dylib
chmod -h og+r libpostscriptlight_mex.dylib
install_name_tool -id libpostscriptlight_mex.dylib libpostscriptlight_mex.5.2.0.dylib
install_name_tool -change ${SLIB}/libz.1.dylib	    ${MLIB}/libz_mex.1.dylib	  libpostscriptlight_mex.5.2.0.dylib
# GMT:
cp ${GLIB}/libgmt.5.2.0.dylib libgmt_mex.5.2.0.dylib
ln -s libgmt_mex.5.2.0.dylib libgmt_mex.5.dylib
chmod -h og+r libgmt_mex.5.dylib
ln -s libgmt_mex.5.dylib     libgmt_mex.dylib
chmod -h og+r libgmt_mex.dylib
install_name_tool -id libgmt_mex.dylib libgmt_mex.5.2.0.dylib
install_name_tool -change ${PLIB}/libnetcdf.7.dylib ${MLIB}/libnetcdf_mex.7.dylib libgmt_mex.5.2.0.dylib
install_name_tool -change ${PLIB}/libgdal.1.dylib   ${MLIB}/libgdal_mex.1.dylib   libgmt_mex.5.2.0.dylib
install_name_tool -change ${PLIB}/libgdal.1.dylib   ${MLIB}/libgdal_mex.1.dylib   libgmt_mex.5.2.0.dylib
install_name_tool -change ${SLIB}/libz.1.dylib	    ${MLIB}/libz_mex.1.dylib	  libgmt_mex.5.2.0.dylib
install_name_tool -change ${GLIB}/libpostscriptlight.5.dylib    ${MLIB}/libpostscriptlight_mex.5.dylib	  libgmt_mex.5.2.0.dylib
cp ${GLIB}/gmt/plugins/supplements.so ${MLIB}/gmt/plugins
cd ${MLIB}/gmt/plugins
install_name_tool -change ${GLIB}/libgmt.5.dylib    ${MLIB}/libgmt_mex.5.dylib    supplements.so
install_name_tool -change ${GLIB}/libpostscriptlight.5.dylib    ${MLIB}/libpostscriptlight_mex.5.dylib    supplements.so
install_name_tool -change ${PLIB}/libnetcdf.7.dylib ${MLIB}/libnetcdf_mex.7.dylib supplements.so
install_name_tool -change ${PLIB}/libgdal.1.dylib   ${MLIB}/libgdal_mex.1.dylib   supplements.so
install_name_tool -change ${SLIB}/libz.1.dylib	    ${MLIB}/libz_mex.1.dylib	  supplements.so
cd ${MLIB}
cd ..
chmod -R og+r mex
cd $HERE
echo "Done"
