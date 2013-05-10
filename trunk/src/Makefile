#	$Id$
#
#	makefile for gmtmex/src directory

include ../config.mk

.SUFFIXES:	.m .$(MEX_EXT)

FLAGS		= $(GMT_INC) $(MEX_INC) -I. -g
ALLLIB     	= $(GMT_LIB) $(MEX_LIB)

PROGS_C		= gmt5.c
PROGS           = $(PROGS_C:.c=.$(MEX_EXT))
PROGS_H		= gmtmex.h
LIB_H		= gmtmex_id.h gmtmex_keys.h gmtmex_progs.h
LIB_C		= gmtmex_parser.c

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
		rm -f $(LIB_H) $(PROGS_H)
clean:	
		rm -f *.o *% core tags $(PROGS)

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

gmt5.$(MEX_EXT):	$(PROGS_C) $(LIB_C) $(LIB_H)
		$(MEX_BLD) $(FLAGS) gmt5.c $(LIB_C) $(ALLLIB)

#-------------------------------------------------------------------------------
#	Build include files
#-------------------------------------------------------------------------------

gmtmex_id.h:	mexproginfo.txt
		echo "enum GMT_prog_enum {" > gmtmex_id.h
		grep -v '^#' mexproginfo.txt | awk '{printf "\tk_%s = %d,\n", $$1, NR-1}' >> gmtmex_id.h
		echo "	k_dummy = -1};" >> gmtmex_id.h

gmtmex_keys.h:	mexproginfo.txt
		echo "static char *keys[] = {" > gmtmex_keys.h
		grep -v '^#' mexproginfo.txt | awk '{printf "\t%s,\n", $$3}' >> gmtmex_keys.h
		echo "	\"\"};" >> gmtmex_keys.h

gmtmex_progs.h:	mexproginfo.txt
		echo "static char *func[] = {" > gmtmex_progs.h
		grep -v '^#' mexproginfo.txt | awk '{printf "\t\"GMT_%s\",\n", $$1}' >> gmtmex_progs.h
		echo "	\"\"};" >> gmtmex_progs.h

