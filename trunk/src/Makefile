#	$Id$
#
#	makefile for gmtmex/src directory

include ../config.mk

INCLUDES	= $(GMT_INC) -I.

PROGS_C		= mexparser.c
PROGS_O         = $(PROGS_C:.c=.o)
PROGS           = $(PROGS_C:.c=)

LIB_C		= gmt_mex.c
LIB_H		= gmt_mex.h
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

clean:	
		rm -f *.o *% *.a core tags $(PROGS)

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

$(PROGS):	$(PROGS_O) $(LIB)
		$(CC) $(LDFLAGS) $@.o $(GMT_LIB) -L../lib -lgmtmex $(LIBS) -o $@

$(LIB_O):	$(LIB_C) $(LIB_H)
