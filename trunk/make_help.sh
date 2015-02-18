#!/bin/bash
#	$Id$
# Build the gmt.m help script from gmt source codes
#
SRC=/Users/pwessel/UH/RESEARCH/CVSPROJECTS/GMTdev/gmt5-dev/trunk/src
cat << EOF > trim.sed
s/{//g
s/}//g
s/,$//g
s/\", /\"\	/g
EOF
V=`gmt-config --version`
Y=`date +%Y`
cat << EOF > /tmp/tmp.m
%
%	GMT - The Generic Mapping Tools, Version $V
%(c) 1991-$Y Paul Wessel, Walter H. F. Smith, R. Scharroo, J. Luis, and F. Wobbe
%
%Supported in part by the US National Science Foundation (http://www.nsf.gov/)
%and volunteers from around the world (see http://gmt.soest.hawaii.edu/).
%
%
%===  GMT core: The main modules of the Generic Mapping Tools  ===
%
%Module name:     Purpose of core module:
% ----------------------------------------------------------------
EOF
grep "\"core\"" $SRC/gmt_core_module.c | grep -v "&" | sed -f trim.sed | awk -F'\t' '{printf "%%%-16s %s\n", $2, $4}' >> /tmp/tmp.m
cat << EOF >> /tmp/tmp.m
%
%===  GMT suppl: The official supplements to the Generic Mapping Tools  ===
%
EOF
for dir in gshhg img meca mgd77 misc potential segy spotter x2sys; do
	cat <<- EOF >> /tmp/tmp.m
	%
	%Module name:     Purpose of $dir module:
	%----------------------------------------------------------------
	EOF
	grep "\"$dir\"" $SRC/gmt_supplements_module.c | sed -f trim.sed | awk -F'\t' '{printf "%%%-16s %s\n", $2, $4}' >> /tmp/tmp.m
done
# Eliminate quotes
sed -e 's/\"//g' /tmp/tmp.m > gmt.m
rm -f trim.sed /tmp/tmp.m