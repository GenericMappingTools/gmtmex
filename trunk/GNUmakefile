#-------------------------------------------------------------------------------
#  $Id: GNUmakefile 389 2015-11-14 06:15:38Z pwessel $
#
#	Guru Makefile for gmt-mex directory
#	
#	!!! THIS MAKEFILE IS FOR GMT-MEX DEVELOPERS ONLY !!!
#
#	This makefile simply makes a tar ball of /opt/gmt.
#
#	Author:	Paul Wessel, SOEST, University of Hawaii
#
#	Date:		15-NOV-2015
#-------------------------------------------------------------------------------
include Makefile

GNUTAR		= $(shell which gnutar || which gtar || which tar)

help::
		@grep '^#!' GNUmakefile | cut -c3-
#!---------------- MAKE HELP FOR GMT-MEX GURUS ----------------
#!
#!opt     : Duplicate active GMT distro to /opt/gmt and re-baptize
#!build   : Configure, build and install the gmt mex files into /opt/gmt/bin
#!tar     : Create a tar ball of the TOOLS source codes
#!update  : Call svn update
#!

opt:
		/Applications/GMT-`gmt --version`.app/Contents/Resources/share/tools/gmt_prepmex.sh

build:
		autoconf
		configure --enable-matlab --enable-debug
		make all
		make install

tar:
		COPYFILE_DISABLE=true $(GNUTAR) --owner 0 --group 0 --mode a=rX,u=rwX \
			-cvjf gmt-`gmt --version`-mex.tar.bz2 /opt/gmt

update:
		svn update
