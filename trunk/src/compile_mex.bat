@echo off
REM ----------------------------------------------------
REM
REM	$Id: compile_mex.bat 113 2013-04-15 21:18:29Z pwessel $
REM
REM
REM	Copyright (c) 1991-2010 by P. Wessel, W. H. F. Smith, R. Scharroo, and J. Luis
REM	See LICENSE.TXT file for copying and redistribution conditions.
REM
REM	This program is free software; you can redistribute it and/or modify
REM	it under the terms of the GNU Lesser General Public License as published by
REM	the Free Software Foundation; version 3 or any later version.
REM
REM	This program is distributed in the hope that it will be useful,
REM	but WITHOUT ANY WARRANTY; without even the implied warranty of
REM	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
REM	GNU Lesser General Public License for more details.
REM
REM	Contact info: gmt.soest.hawaii.edu
REM --------------------------------------------------------------------
REM --------------------------------------------------------------------------------------
REM
REM This is a compile batch that builds all GMT MEXs. Contrary to the 'mex' command it doesn't
REM need you to setup a compiler within MATLAB, which means you can use any of the MS or
REM Intel compilers (a luxury that you don't have with the 'mex' command).
REM
REM If a WIN64 version is targeted than both GMT & netCDF Libs must have been build in 64-bits as well.
REM
REM
REM Usage: open the command window set up by the compiler of interest (were all vars are already set)
REM	   and run this batch from there.
REM 	   NOTE: you must make some edits to the setup below.
REM
REM --------------------------------------------------------------------------------------

REM ------------- Set the compiler (set to 'icl' to use the Intel compiler) --------------
SET CC=cl
REM --------------------------------------------------------------------------------------

REM If set to "yes", linkage is done againsts ML6.5 Libs
SET R13="no"

REM Set it to 32 or 64 to build under 64-bits or 32-bits respectively.
SET BITS=64

IF %R13%=="yes" SET BITS=32

REM Options are "dll", "mexw32" (recent ML version scream when they find .dll) or "mexw64" (when BITS=64)
SET MEX_EXT="mexw32"

REM
REM Set to "yes" if you want to build a debug version
SET DEBUG="yes"
REM
SET LDEBUG=
IF %DEBUG%=="yes" SET LDEBUG=/debug

REM ------------------ Sets the MATLAB libs and include path ----------------------------
IF %R13%=="yes" (

SET MATLIB=C:\SVN\pracompila\MAT65\lib\win32\microsoft
SET MATINC=C:\SVN\pracompila\MAT65\include
SET _MX_COMPAT=
SET MEX_EXT="dll"

) ELSE (

IF %BITS%==64 (
SET MATLIB=C:\SVN\pracompila\ML2010a_w64\lib\win64\microsoft 
SET MATINC=C:\SVN\pracompila\ML2010a_w64\include
SET _MX_COMPAT=-DMX_COMPAT_32
SET MEX_EXT="mexw64"

) ELSE (

SET MATLIB=C:\SVN\pracompila\ML2009b_w32\lib\win32\microsoft 
SET MATINC=C:\SVN\pracompila\ML2009b_w32\include
SET _MX_COMPAT=-DMX_COMPAT_32
SET MEX_EXT="mexw32"
) )

REM -------------- Set GMT & NetCDF lib and include ----------------------------
IF %BITS%==64 (

SET  NETCDF_LIB=C:\programs\compa_libs\netcdf\compileds\VC10_64\lib\netcdf.lib
SET     GMT_LIB=c:\progs_cygw\GMTdev\gmt5\trunk\WIN%BITS%\lib\gmt.lib
SET    GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_64\lib\gdal_i.lib
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" amd64

) ELSE (

SET  NETCDF_LIB=C:\programs\compa_libs\netcdf\compileds\VC10_32\lib\netcdf.lib
SET     GMT_LIB=c:\progs_cygw\GMTdev\gmt5\trunk\WIN%BITS%\lib\gmt.lib
SET    GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_32\lib\gdal_i.lib
call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\vcvarsall.bat" x86

)

SET  NETCDF_INC=C:\programs\compa_libs\netcdf\compileds\VC10_32\include
SET     GMT_INC=c:\progs_cygw\GMTdev\gmt5\trunk\WIN%BITS%\include
SET    GMT_INC2=c:\progs_cygw\GMTdev\gmt5\trunk\src
SET    GDAL_INC=c:\programs\GDALtrunk\gdal\compileds\VC10_32\include
REM ----------------------------------------------------------------------------

REM ____________________________________________________________________________
REM ___________________ STOP EDITING HERE ______________________________________


SET COMPFLAGS=/Zp8 /GR /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD
SET OPTIMFLAGS=/Ox /Oy- /DNDEBUG
IF %DEBUG%=="yes" SET OPTIMFLAGS=/Z7

IF %BITS%==64 SET arc=X64
IF %BITS%==32 SET arc=X86
SET LINKFLAGS=/dll /export:mexFunction /LIBPATH:%MATLIB% libmx.lib libmex.lib libmat.lib /MACHINE:%arc% kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /incremental:NO %LDEBUG%


REM -------------- GENERATE THESE GUYS --------------------------------------------------------------------
REM echo "enum GMT_prog_enum {" > gmtmex_id.h
REM grep -v "#" mexproginfo.txt | gawk "{printf \"\tk_%s = %d,\n\", $1, NR-1}" >> gmtmex_id.h
REM echo "k_dummy = -1};" >> gmtmex_id.h

REM echo "static char *keys[] = {" > gmtmex_keys.h
REM grep -v "#" mexproginfo.txt | gawk "{printf \"\t%s,\n\", $$3}" >> gmtmex_keys.h
REM echo "\"\"};" >> gmtmex_keys.h

REM -------------------------------------------------------------------------------------------------------


REM -------------------------------------------------------------------------------------------------------
%CC% /c -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% -I%GMT_INC% %OPTIMFLAGS% %_MX_COMPAT% -DLIBRARY_EXPORTS -DGMT_MATLAB gmtmex_parser.c gmt5.c gmt.c


rem link  /out:"mexparser.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% /implib:templib.x gmtmex_parser.obj gmt5.obj
rem %CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% -I%GMT_INC% %OPTIMFLAGS% %_MX_COMPAT% -DLIBRARY_EXPORTS -DGMT_MATLAB /DFUNC=GMT_grdproject /Fegrdproject mexprogram.c

link  /out:"gmt.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% /implib:templib.x gmtmex_parser.obj gmt.obj
REM -------------------------------------------------------------------------------------------------------


del *.obj *.exp templib.x

pause
