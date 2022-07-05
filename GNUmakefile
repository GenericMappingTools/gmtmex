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
#	Date:		1-NOV-2021
#-------------------------------------------------------------------------------
include Makefile

GNUTAR		= $(shell which gnutar || which gtar || which tar)

help::
		@grep '^#!' GNUmakefile | cut -c3-
#!---------------- MAKE HELP FOR GMT-MEX GURUS ----------------
#!
#!opt     : Duplicate active GMT distro to /opt/gmt and re-baptize
#!build   : Configure, build and install the gmt mex files into /opt/gmt/bin
#!tar     : Create a tar ball of the mex/gmt biinary distro
#!update  : Call git pull
#!wipe    : Remove mex-* tarballs
#!latest-config : Update the configure include files
#!

opt:
		@echo "[Running `ls /Applications/GMT-6.4.?.app/Contents/Resources/share/tools/gmt_prepmex.sh | tail -1`]"; echo ""
		@`/Applications/GMT-6.4.?.app/Contents/Resources/share/tools/gmt_prepmex.sh | tail -1`

latest-config:
		curl "http://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub;hb=HEAD" -s -R -o config.sub
		curl "http://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess;hb=HEAD" -s -R -o config.guess

build:
		autoconf
		gmtswitch /opt/gmt
		configure --enable-matlab --enable-debug
		make all
		make install

tar:
		COPYFILE_DISABLE=true $(GNUTAR) --owner 0 --group 0 --mode a=rX,u=rwX --absolute-names \
			-cvjf mex-gmt-`gmt --version`-darwin-x84_64.tbz /opt/gmt

update:
		git pull

wipe:
		rm -f mex-gmt-*-darwin-x84_64.tbz
