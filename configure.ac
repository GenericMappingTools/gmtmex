dnl -----------------------------------------------------------------
dnl
dnl     This program is free software; you can redistribute it and/or modify
dnl     it under the terms of the GNU General Public License as published by
dnl     the Free Software Foundation; version 2 of the License.
dnl
dnl     This program is distributed in the hope that it will be useful,
dnl     but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl     GNU General Public License for more details.
dnl
dnl     Contact info: Paul Wessel, SOEST/UHM, pwessel@hawaii.edu
dnl ----------------------------------------------------------------------------
dnl
dnl Run GNU autoconf to produce a configure script from this input.
dnl -----------------------------------------------------------------
dnl
dnl
dnl -----------------------------------------------------------------
dnl Since gmtmex requires GMT5 libs we insist that the GMT5
dnl executables must be in the current path; otherwise abort.
dnl 
dnl We also need to gather the following information
dnl 1. Location of GMT5 include and library files
dnl -----------------------------------------------------------------
dnl
AC_INIT([gmt-mex],[6],[pwessel@hawaii.edu])
dnl
dnl -----------------------------------------------------------------
dnl The remainder should need no editing
dnl -----------------------------------------------------------------
dnl
AC_CONFIG_SRCDIR(src/gmtmex.h)
AC_CONFIG_AUX_DIR(`pwd`)
AC_PREFIX_PROGRAM(gmt)	
m4_warn([obsolete],
[AC_CYGWIN is obsolete: use AC_CANONICAL_HOST and check if $host_os
matches *cygwin*])dnl
AC_CANONICAL_HOST
case $host_os in
  *cygwin* ) CYGWIN=yes;;
	 * ) CYGWIN=no;;
esac

AC_CANONICAL_TARGET
AC_LANG([C])
AC_PROG_CC
AC_PROG_CPP
AC_PREFIX_DEFAULT(`pwd`)
AC_PATH_XTRA
dnl
dnl -----------------------------------------------------------------
dnl gmtmex-specific make macros
dnl -----------------------------------------------------------------
dnl
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LD_OPT)
AC_SUBST(LDFLAGS)
AC_SUBST(GMT_LIB)
AC_SUBST(GMT_INC)
AC_SUBST(GMT_BIN)
AC_SUBST(DEP_LIBS)
AC_SUBST(LIBS)
AC_SUBST(LIBEXT)
AC_SUBST(SL)
AC_SUBST(SL_VERSION)
AC_SUBST(prefix)
AC_SUBST(rootdir)
AC_SUBST(MEX_BLD)
AC_SUBST(MEX_LIB)
AC_SUBST(MEX_INC)
AC_SUBST(MEX_EXT)
AC_SUBST(MEX_OUT)
AC_SUBST(MEX_XDIR)
AC_SUBST(MEX_MDIR)
AC_SUBST(GS_PATH)
dnl
dnl -----------------------------------------------------------------
dnl Special configure options for gmtmex installation
dnl -----------------------------------------------------------------
dnl
AC_ARG_WITH(gmt_config, [  --with-gmt-config=      Full path to gmt-config [[auto]]])
AC_ARG_ENABLE(matlab,   [  --enable-matlab=dir     Build GMT-matlab interface, give directory with MATLAB distribution [[auto]]])
AC_ARG_ENABLE(octave,   [  --enable-octave=dir     Build GMT-octave interface, give directory with Octave distribution [[auto]] (disables MATLAB interface)])
AC_ARG_ENABLE(octmode,  [  --enable-octmode=mode   Append oct or mex [[auto]]])
AC_ARG_ENABLE(mex_mdir, [  --enable-mex-mdir=dir   Directory in which to install GMT-octave/MATLAB .m files [[auto]]])
AC_ARG_ENABLE(mex_xdir, [  --enable-mex-xdir=dir   Directory in which to install GMT-octave/MATLAB mex files [[auto]]])
AC_ARG_ENABLE(debug,    [  --enable-debug          Compile for debugging instead of optimizing code])
AC_ARG_ENABLE(rpath,	[  --disable-rpath         Do not hardcode runtime library paths])
AC_ARG_ENABLE(shared,   [  --enable-shared         Build shared (dynamic) libraries instead of static])
AC_ARG_VAR(GMT_INC,Location of GMT5 headers (compile-time))
AC_ARG_VAR(GMT_LIB,Location of GMT5 library (compile-time))
dnl
dnl -----------------------------------------------------------------
dnl Check for ar, csh/bash, ld, awk, and ranlib
dnl -----------------------------------------------------------------
dnl
AC_CHECK_PROG(GMT_in_path, gmt, 1, 0)
if test "X$GMT_in_path" = "X0" ; then
	AC_MSG_RESULT(GMT5 executables are not in your path - aborting)
	exit 1
else
	gmajor=`gmt --version | awk '{print substr($1,1,1)}'`
	gminor=`gmt --version | awk '{print substr($1,3,1)}'`
	if test "X$gmajor" = "X6"; then
		AC_MSG_RESULT(GMT6.x.x or higher executables found)
	elif test "X$gmajor" = "X5" && test $gminor -ge 3; then
		AC_MSG_RESULT(GMT5.3.x or higher executables found)
	else
		AC_MSG_RESULT(GMT5.3.x executables are not in your path - aborting)
		exit 1
	fi
fi
AC_CHECK_PROGS(PAGER, less more, cat)
AC_CHECK_PROG(AR, ar, ar, :)
AC_CHECK_PROG(CSH, csh, csh, bash, )
AC_CHECK_PROG(LD, ld, ld, :)
AC_CHECK_PROG(AS, as, as, :)
AC_PROG_LN_S
AC_PROG_INSTALL
AC_PROG_AWK
AC_PROG_RANLIB
AC_CHECK_SIZEOF(long)
dnl
dnl -----------------------------------------------------------------
dnl Find full path to gs
dnl -----------------------------------------------------------------
AC_PATH_PROG([GS_EXE], [gs], [$gs_exe])
if test -x "$GS_EXE" ; then
	GS_PATH=$GS_EXE
fi	
dnl
dnl -----------------------------------------------------------------
dnl Get current OS and chip type
dnl -----------------------------------------------------------------
dnl
os=`uname -s`
if test "X$os" = "XDarwin" ; then
	osr=`uname -r | cut -d. -f1`
fi
chip=`uname -p`
mchip=`uname -m`
Sun_rev=`uname -r | grep '^5.'`
rootdir=`pwd`
dnl
dnl -----------------------------------------------------------------
dnl Check for libraries.  Most installs will need to have
dnl -lm but under OS X it is not needed.  Thus, we look
dnl for sin in -lm which should set -lm on most platforms.
dnl We also look for sincos in -m and -lsunmath for Solaris.
dnl -----------------------------------------------------------------
dnl
AC_SEARCH_LIBS(sin, m)
dnl
dnl -----------------------------------------------------------------
dnl We use gmt-config to set GMT5 paths and settings
dnl -----------------------------------------------------------------
dnl
if test "X$with_gmt_config" != "Xno" ; then
	AC_PATH_PROG([GMT_CONF], [gmt-config], [$with_gmt_config])
fi
if test -x "$GMT_CONF" && test "X$GMT_INC" = "X" && test "X$GMT_LIB" = "X" ; then
	AC_MSG_CHECKING(for GMT includes)
	GMT_INC=`$GMT_CONF --cflags`
	AC_MSG_RESULT($GMT_INC)
	AC_MSG_CHECKING(for GMT libs)
	GMT_LIB=`$GMT_CONF --libs`
	DEP_LIBS=`$GMT_CONF --dep-libs`
	AC_MSG_RESULT($GMT_LIB)
	AC_MSG_CHECKING(for GMT version)
	GMT_VER=`$GMT_CONF --version`
	AC_MSG_RESULT($GMT_VER)
	AC_MSG_CHECKING(if GMT was built with GDAL)
	GMT_GDAL=`$GMT_CONF --has-gdal`
	AC_MSG_RESULT($GMT_GDAL)
	AC_MSG_CHECKING(for GMT bin directory)
	GMT_BIN=`gmt --show-bindir`
	AC_MSG_RESULT($GMT_BIN)
	GMT_LIB_PATH=`echo $GMT_LIB | sed 's/^-L//;s/\ .*//'`
	rpath="$GMT_LIB_PATH"
	AC_MSG_CHECKING(for 32 vs 64-bit GMT5 installation)
	bits=`$GMT_CONF --bits`
	GCC_64="-m${bits}"
	LDFLAGS="$LDFLAGS -m${bits}"
	AC_MSG_RESULT(found ${bits}-bit)
fi
dnl
dnl -----------------------------------------------------------------
dnl Allow user to set --enable-debug to pass -g flag
dnl -----------------------------------------------------------------
dnl
AC_MSG_CHECKING(if we should compile with debugging switch -g)
CFLAGS=$ac_save_CFLAGS
if test "X$enable_debug" = "Xyes" ; then
	AC_MSG_RESULT(yes)
	CFLAGS="$CFLAGS -g"
	MFLAGS="-g"
else
	CFLAGS=`echo $CFLAGS | sed 's/-g//'`
	if test "$os" = "Darwin" ; then
		LDFLAGS="$LDFLAGS -Wl,-x"
	else
		LDFLAGS="$LDFLAGS -s"
	fi
	if test "$GCC" = "yes" || test "$os" = "Linux" || test "$os" = "AIX"; then
		CFLAGS="$CFLAGS -O2"
	elif test "$chip" = "alpha" || test "$os" = "AUX" || test "$os" = "HP-UX"; then
		CFLAGS="$CFLAGS"
	else
		CFLAGS="$CFLAGS -O"
	fi
	AC_MSG_RESULT(no)
fi
dnl
dnl -----------------------------------------------------------------
dnl Determine LD_OPT for this platform
dnl -----------------------------------------------------------------
dnl
shared_unknown=no
LD_OPT=
LD=ld
SL=so
LIBEXT='$(SL)'
SL_VERSION='$(SL)'
if test "X$enable_shared" = "Xyes" ; then
	AC_MSG_CHECKING(how to make a shared object)
	if test "X$CYGWIN" = "Xyes" ; then
		SL=dll.a
		SL_VERSION=$LIB_MAJOR_VERSION.$SL
		LD='$(CC)'
		LD_OPT="-shared"
	elif test "$os" = "SunOS" ; then
		if test "$Sun_rev" = "" ; then	# SunOS
			LD_OPT="-Bdynamic"
		else				# Solaris
			LD_OPT="-G"
		fi
	elif test "$os" = "Linux" || test "$os" = "GNU/kFreeBSD" ; then	# Need optimization when doing shared
		LD='$(CC)'
		SL_VERSION=$SL.$LIB_MAJOR_VERSION
		LD_OPT='-shared -Wl,-soname=$*.'$SL_VERSION
	elif test "$os" = "FreeBSD" || test "$os" = "NetBSD" ; then
		SL_VERSION=$SL.$LIB_MAJOR_VERSION
		LD_OPT='-shared -Wl,-soname=$*.'$SL_VERSION
	elif test "$os" = "OSF1" ; then
		LD_OPT='-shared -expect_unresolved "*"'
	elif test "$os" = "IRIX64" ; then
		LD_OPT="-shared"
		# Need to remove -L/usr/lib to avoid mixing o32 and n32
		X_LIBS=
	elif test "$os" = "Darwin" ; then	# MacOS X
		SL=dylib
		SL_VERSION=$LIB_MAJOR_VERSION.$SL
		LD='$(CC)'
		LD_OPT='-install_name ${prefix}/lib/$*.'$SL_VERSION' -compatibility_version 4.0.0 -current_version '$LIB_VERSION' -dynamiclib '$OSX_L
		if test $osr -eq 8 ; then	# Tiger needs special option for shared libs
			CFLAGS="$CFLAGS -fno-common"
		fi
	elif test "$os" = "HP-UX" ; then
		LD_OPT="-b"
		SL=sl
	else	# Not yet implemented
		AC_MSG_RESULT(not supported for gmtmex under $os - building static libs instead)
		enable_shared=no
		shared_unknown=yes
	fi
fi
if test "X$enable_shared" = "Xyes" ; then
	rpath="$rpath"':$(libdir)'
	AC_MSG_RESULT($LD $LD_OPT)
	AC_MSG_CHECKING(the file extension for shared libraries)
	AC_MSG_RESULT($SL)
else
	LIBEXT='a'
fi
dnl
dnl -----------------------------------------------------------------
dnl Determine LDFLAGS for this platform
dnl -----------------------------------------------------------------
dnl
AC_MSG_CHECKING(linker options)
flag=
if test "X$enable_rpath" = "Xno" ; then
	flag=
elif test "$os" = "SunOS" ; then
	if test "$Sun_rev" = "" ; then
		flag="-L"
	else
		flag="-R,"
	fi
elif test "$os" = "Linux" || test "$os" = "FreeBSD" || test "$os" = "NetBSD" || test "$os" = "OSF1" || test "$os" = "IRIX64" ; then
	flag="-rpath,"
elif test "$os" = "HP-UX" ; then
	flag="+b,"
fi
if test "X$flag" != "X" ; then
	LDFLAGS="$LDFLAGS -Wl,$flag$rpath"
fi
AC_MSG_RESULT($LDFLAGS)
dnl
dnl -----------------------------------------------------------------
dnl Determine CFLAGS for this platform
dnl -----------------------------------------------------------------
dnl
AC_MSG_CHECKING(compiler options for GMT)
if test "$GCC" = "yes"; then		# GNU cc options
	if test "$chip" = "alpha" ; then	# Alpha PC
		CFLAGS="$CFLAGS -mieee -Wall $GCC_64"
	else				# MacOSX, Linux, ...
		CFLAGS="$CFLAGS -Wall $GCC_64"
	fi
	CFLAGS="$CFLAGS -fPIC -fno-strict-aliasing -std=c99"
elif test "$os" = "SunOS" ; then	# Sun Workshop CC options
	CFLAGS="$CFLAGS -Xc"
elif test "$os" = "HP-UX" ; then	# HP cc options
	CFLAGS="$CFLAGS -Aa +Ofltacc +O1 +Olibcalls"
	if test "$enable_shared" = "yes" ; then
		CFLAGS="$CFLAGS +z"
	fi
elif test "$os" = "OSF1" ; then		# Dec Alpha OSF cc options
	CFLAGS="$CFLAGS -Olimit 1500 -ieee_with_inexact"
elif test "$os" = "ULTRIX" ; then	# Dec Ultrix cc options
	CFLAGS="$CFLAGS -Olimit 1500"
fi
AC_MSG_RESULT($CFLAGS)
dnl
dnl -----------------------------------------------------------------
dnl Determine CPPFLAGS for this platform
dnl -----------------------------------------------------------------
dnl
AC_MSG_CHECKING(compiler preprocessor flags for gmtmex)
if test "$os" = "AUX" ; then		# Apple A/UX cc options
	CPPFLAGS="$CPPFLAGS -Daux"
fi
AC_MSG_RESULT($CPPFLAGS)
dnl
dnl -----------------------------------------------------------------
dnl Check if we have GNU make. If not, add CPPFLAGS to CFLAGS
dnl -----------------------------------------------------------------
dnl
AC_MSG_CHECKING(if CPPFLAGS need to be added to CFLAGS)
${MAKE:-make} -p -f /dev/null 2>&1 | grep CPPFLAGS > /dev/null
if test "$?" = "0" ; then
	AC_MSG_RESULT(no)
else
	CFLAGS="$CFLAGS "'$(CPPFLAGS)'
	AC_MSG_RESULT(yes)
fi
dnl
dnl
dnl -----------------------------------------------------------------
dnl Determine Octave/MATLAB mex extensions and dirs for this platform
dnl -----------------------------------------------------------------
dnl
if test "X$MATLAB" = "XNONE" ; then
	AC_MSG_RESULT([Octave/MATLAB mex supplement will be skipped])
elif test ! "X$enable_matlab" = "X"; then
	AC_MSG_CHECKING(Compiler options for MATLAB)
	MEX_BLD="$MATLAB/bin/mex -DGMT_MATLAB $MFLAGS"
	MEX_OUT='-output'
	CFLAGS="$CFLAGS -DGMT_MATLAB"
	if test "$os" = "Linux" ; then		# Linux systems
		if test "$ac_cv_sizeof_long" = "8" ; then
			MEX_EXT='mexa64'
			MATLAB_MEX='glnxa64'
		else
			MEX_EXT='mexglx'
			MATLAB_MEX='glnx86'
		fi
	elif test "$os" = "SunOS" ; then	# Sun Workshop
		MEX_EXT='mexsol'
		MATLAB_MEX='sol2'
	elif test "$os" = "IRIX64" ; then	# SGI Irix 64
		MEX_EXT='mexsg64'
	elif test "$os" = "IRIX" ; then		# SGI Irix
		MEX_EXT='mexsg'
	elif test "$os" = "HP-UX" ; then	# HP
		MEX_EXT='mexhpux'
	elif test "$os" = "AIX" ; then		# IBM AIX
		MEX_EXT='mexrs6'
		MATLAB_MEX='ibm_rs'
	elif test "$os" = "OSF1" ; then		# Dec Alpha
		MEX_EXT='mexds'
	elif test "$os" = "ULTRIX" ; then	# Dec Ultrix
		MEX_EXT='mexvx'
	elif test "$os" = "Darwin" ; then	# Mac OS X
		dir=`ls -t /Applications | grep MATLAB | head -n 1`
		if test "X$dir" = "X" ; then
			MATLAB=NONE
		else
			MATLAB=/Applications/$dir
		fi
		MEX_BLD="xcrun clang -undefined error -arch ${mchip} -bundle -DGMT_MATLAB $MFLAGS"
		MEX_OUT='-o'
		if test "$chip" = "i386" ; then		#Intel Mac
			if test "$ac_cv_sizeof_long" = "8" ; then
				MEX_EXT='mexmaci64'
				MATLAB_MEX='maci64'
			else
				MEX_EXT='maci64'
				MATLAB_MEX='maci'
			fi
		elif test "$chip" = "arm" ; then		#Apple Silicon Mac 64-bit
			MEX_EXT='mexmaca64'
			MATLAB_MEX='maca64'
		else	#Legacy PPC
			MEX_EXT='mexmac'
			MATLAB_MEX='mac'
		fi
	else					# Other
		MEX_EXT='mex'
	fi
	if test "X$MATLAB_MEX" = "X" ; then
		# Users: pls supply your input
		AC_MSG_RESULT([MATLAB mex-file extension not known - guessed to be mex])
		AC_MSG_RESULT([If incorrect, please modify makefile manually and post an issue on the GMT/MEX wiki])
		MATLAB_MEX='mex'
	fi
	MEX_LIB="-L$MATLAB/bin/$MATLAB_MEX -lmx -lmex"
	MEX_INC="-I$MATLAB/extern/include"
	if test "X$enable_mex_mdir" = "X" || test "X$enable_mdir" = "Xyes" ; then
		MEX_MDIR=$MATLAB/toolbox/local
	else
		MEX_MDIR=$enable_mex_mdir
	fi
	if test "X$enable_mex_xdir" = "X" || test "X$enable_xdir" = "Xyes" ; then
		MEX_XDIR=$MATLAB/toolbox/local
	else
		MEX_XDIR=$enable_mex_xdir
	fi
	AC_MSG_RESULT([$MEX_BLD $MEX_LIB, extension .$MEX_EXT])
	AC_MSG_CHECKING(MATLAB install dir for .m files)
	AC_MSG_RESULT([$MEX_MDIR])
	AC_MSG_CHECKING(MATLAB install dir for .mex files)
	AC_MSG_RESULT([$MEX_XDIR])
else
	if test "X$enable_octmode" = "Xyes" || test "X$enable_octmode" = "Xoct" ; then
		MEX_EXT='oct'
		MEX_ARG=
		OCT_ARG=GMT_OCTOCT
	else
		MEX_EXT='mex'
		MEX_ARG='-mex'
		OCT_ARG=GMT_OCTMEX
	fi
	AC_MSG_CHECKING(Compiler options for Octave)
	CFLAGS="$CFLAGS -D$OCT_ARG"
	MEX_BLD="mkoctfile -D$OCT_ARG $MEX_ARG $MFLAGS"
	MEX_LIB=
	MEX_INC=-I`octave-config --print OCTINCLUDEDIR`
	if test "X$enable_mex_mdir" = "X" || test "X$enable_mdir" = "Xyes" ; then
		MEX_MDIR=`octave-config --m-site-dir`
	else
		MEX_MDIR=$enable_mex_mdir
	fi
	if test "X$enable_mex_xdir" = "X" || test "X$enable_xdir" = "Xyes" ; then
		MEX_XDIR=`octave-config --oct-site-dir`
	else
		MEX_XDIR=$enable_mex_xdir
	fi
	AC_MSG_RESULT([$MEX_BLD])
	AC_MSG_CHECKING(Octave extension)
	AC_MSG_RESULT([.$MEX_EXT])
	AC_MSG_CHECKING(Octave include dir)
	AC_MSG_RESULT([$MEX_INC])
	AC_MSG_CHECKING(Octave system install dir for .m files)
	AC_MSG_RESULT([$MEX_MDIR])
	AC_MSG_CHECKING(Octave system install dir for .$MEX_EXT files)
	AC_MSG_RESULT([$MEX_XDIR])
fi
if test "X$GMT_GDAL" = "yes" ; then
	CFLAGS="$CFLAGS -DHAVE_GDAL"
fi
dnl -----------------------------------------------------------------
dnl Update the Make macros
dnl -----------------------------------------------------------------
dnl
AC_CONFIG_FILES([config.mk])
AC_OUTPUT
dnl
