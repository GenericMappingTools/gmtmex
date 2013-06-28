#	$Id$
#
#	makefile for gmtmex/src directory

include ../config.mk

.SUFFIXES:	.m .$(MEX_EXT)

FLAGS		= $(GMT_INC) $(MEX_INC) -I.
ALLLIB     	= $(GMT_LIB) $(MEX_LIB)

PROGS_C		= gmt.c
PROGS           = $(PROGS_C:.c=.$(MEX_EXT))
PROGS_H		= gmtmex.h
LIB_H		= gmtmex_keys.h gmtmex_modules.h
LIB_C		= gmtmex_parser.c
LIB_O		= $(LIB_C:.c=.o)

#-------------------------------------------------------------------------------
#	software targets
#-------------------------------------------------------------------------------

all:		$(PROGS)

install:	all
		$(INSTALL) -d $(bindir)
		$(INSTALL) $(PROGS) $(bindir)

uninstall:
		cd $(bindir); rm -f $(PROGS)

spotless::	clean
		rm -f $(LIB_H)
clean:	
		rm -f *.o *% core tags $(PROGS) gmt_mextest

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

gmt5.$(MEX_EXT):	$(PROGS_C) $(LIB_C) $(LIB_H)
		$(MEX_BLD) $(FLAGS) gmt5.c $(LIB_C) $(ALLLIB)

gmt.$(MEX_EXT):	$(PROGS_C) $(LIB_C) $(LIB_H)
		$(MEX_BLD) $(FLAGS) gmt.c $(LIB_C) $(ALLLIB)

gmt_mextest:	gmt_mextest.o $(LIB_O)
		$(CC) $(LDFLAGS) -DNO_MEX gmt_mextest.o $(LIB_O) $(GMT_LIB) $(LIBS) -o $@

$(LIB_O):	$(LIB_C) $(LIB_H)
		$(CC) $(CFLAGS) -DNO_MEX $(LIB_C) -c -o $@

gmt_mextest.o:	gmt_mextest.c
		$(CC) $(CFLAGS) -DNO_MEX gmt_mextest.c -c -o $@
	
#-------------------------------------------------------------------------------
#	Build include files
#-------------------------------------------------------------------------------

gmtmex_keys.h:	mexproginfo.txt
		echo "static char *keys[] = {" > gmtmex_keys.h
		grep -v '^#' mexproginfo.txt | sort -k1 | awk '{printf "\t%s,\n", $$3}' >> gmtmex_keys.h
		echo "	\"\"};" >> gmtmex_keys.h

gmtmex_modules.h:	mexproginfo.txt
		echo "static char *module_name[] = {" > gmtmex_modules.h
		grep -v '^#' mexproginfo.txt | sort -k1 | awk '{printf "\t\"%s\",\n", $$1}' >> gmtmex_modules.h
		echo "	\"\"};" >> gmtmex_modules.h
		grep -v '^#' mexproginfo.txt | sort -k1 | wc -l | awk '{printf "\n#define N_GMT_MODULES %d\n", $$1}' >> gmtmex_modules.h

