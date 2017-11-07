#       $Id$
#
#       makefile for top gmtmex directory

sinclude config.mk

all:
	cd src; $(MAKE) all

info:
	bash make_help.sh

test:
	cd src; $(MAKE) gmt_mextest

install:
	cd src; $(MAKE) install

uninstall:
	cd src; $(MAKE) uninstall

clean:
	cd src; $(MAKE) clean

spotless::
	cd src; $(MAKE) spotless
	rm -rf config.log config.status config.mk configure autom4te.cache
