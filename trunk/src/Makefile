#	$Id$
#
#	makefile for gmtmex/src directory

include ../config.mk

INCLUDES	= $(GMT_INC) $(MEX_INC) -I.

PROGS_C		= gmt5.c
PROGS_O         = $(PROGS_C:.c=.o)
PROGS           = $(PROGS_C:.c=)
PROGS_H		= gmtmex.h
LIB_H		= gmtmex_id.h gmtmex_keys.h gmtmex_progs.h
LIB_C		= gmtmex_parser.c
LIB_O		= $(LIB_C:.c=.o)
LIB		= libgmtmex.$(LIBEXT)

#-------------------------------------------------------------------------------
#	software targets
#-------------------------------------------------------------------------------

lib:		$(LIB)

all:		$(PROGS) lib

libgmtmex.$(LIBEXT):	$(LIB_O)
		$(AR) cvur $@ $?
		$(RANLIB) $@

install:	all
		$(INSTALL) -d $(bindir)
		$(INSTALL) $(PROGS) $(bindir)
		$(INSTALL) -d ../lib
		$(INSTALL) -m 644 $(LIB:.$(LIBEXT)=.a) ../lib
		if [ ! $(LIBEXT) = a ]; then $(INSTALL) -m 755 $(LIB) ../lib; fi

uninstall:
		cd $(bindir); rm -f $(PROGS)

spotless::	clean
		rm -f $(LIB_H) $(PROGS_H)
clean:	
		rm -f *.o *% *.a core tags $(PROGS)

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

$(PROGS):	$(PROGS_O) $(LIB)
		$(CC) $(LDFLAGS) $@.o $(GMT_LIB) -L. -lgmtmex $(LIBS) -o $@

$(LIB_O):	$(LIB_C) $(LIB_H)

$(PROGS_O):	$(PROGS_H)

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

