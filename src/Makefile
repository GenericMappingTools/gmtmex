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

clean:	
		rm -f *.o *% core tags $(PROGS) gmt_mextest

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

gmt5.$(MEX_EXT):	$(PROGS_C) $(LIB_C)
		$(MEX_BLD) $(FLAGS) gmt5.c $(LIB_C) $(ALLLIB)

gmt.$(MEX_EXT):	$(PROGS_C) $(LIB_C)
		$(MEX_BLD) $(FLAGS) gmt.c $(LIB_C) $(ALLLIB)

gmt_mextest:	gmt_mextest.o $(LIB_O)
		$(CC) $(LDFLAGS) -DNO_MEX gmt_mextest.o $(LIB_O) $(GMT_LIB) $(LIBS) -o $@

$(LIB_O):	$(LIB_C)
		$(CC) $(CFLAGS)  $(GMT_INC) -DNO_MEX $(LIB_C) -c -o $@

gmt_mextest.o:	gmt_mextest.c
		$(CC) $(CFLAGS)  $(GMT_INC) -DNO_MEX gmt_mextest.c -c -o $@
