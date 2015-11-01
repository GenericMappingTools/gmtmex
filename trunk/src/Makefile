#	$Id$
#
#	makefile for gmtmex/src directory

include ../config.mk

.SUFFIXES:	.m .$(MEX_EXT)

FLAGS		= $(GMT_INC) $(MEX_INC)
ALLLIB     	= $(GMT_LIB) $(MEX_LIB)

PROGS_C		= gmtmex.c
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
		rm -f *.o *% core tags gmtmex.mex* gmtmex.oct gmt_opttest

#-------------------------------------------------------------------------------
#	program rules
#-------------------------------------------------------------------------------

gmt5.$(MEX_EXT):	$(PROGS_C) $(LIB_C)
		$(MEX_BLD) $(FLAGS) gmt5.c $(LIB_C) $(ALLLIB)

gmtmex.$(MEX_EXT):	$(PROGS_C) $(LIB_C)
		$(MEX_BLD) $(FLAGS) gmtmex.c $(LIB_C) $(ALLLIB)

gmt_opttest:	gmt_opttest.o
		$(CC) $(LDFLAGS) gmt_opttest.o $(GMT_LIB) $(LIBS) -o $@

$(LIB_O):	$(LIB_C)
		$(CC) $(FLAGS) $(CFLAGS) $(GMT_INC) -DNO_MEX $(LIB_C) -c -o $@

gmt_opttest.o:	gmt_opttest.c
		$(CC) $(CFLAGS) $(GMT_INC) gmt_opttest.c -c -o $@
