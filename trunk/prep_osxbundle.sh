#!/bin/bash
# $Id: prep_osx.sh 312 2015-04-11 03:38:12Z pwessel $
# Rebaptizing of shared libraries for OSX as distributed by
# the GMT bundle.  We duplicate and rebaptize all the libs
# so gmtmex can link with the libXgmt library instead.
# We also duplicate gs and its library and rebaptize the
# two libraries it uses that would clash with ML versions.
# Run from /Applications/GMT-x.x.x/Contents/Resources/lib
#-----------------------------------------------------------------
# Make shadow directory mex that will contain rebaptized libs:
mkdir -p mex/gmt/plugins
# Find a list of all libs shipped with the OSX bundle:
ls *.dylib | egrep -v 'libgmt.dylib|libpsl.dylib' > /tmp/l.lis
# For each, duplicate into mex but add a leading X to each name
while read lib; do
	new=`echo $lib | awk '{printf "libX%s\n", substr($1,4)}'`
	cp $lib mex/$new
done < /tmp/l.lis
# Same for the supplement shared plugin
cp gmt/plugins/supplements.so mex/gmt/plugins/supplements.so
cd mex
ls *.dylib > /tmp/l.lis
# For all libs in mex, change internal references to contain the leading "X"
while read lib; do
	echo
	echo "--> $lib"
	echo
	otool -L $lib | grep executable_path | awk '{print $1}' > /tmp/t.lis
	let k=1
	while read old; do
		new=`echo $old | awk -F/ '{printf "libX%s\n", substr($NF,4)}'`
		if [ $k -eq 1 ]; then # Do the id change
			was=`echo $lib | awk -F/ '{print substr($1,4)}'`
			install_name_tool -id @rpath/$new $lib
		else
			install_name_tool -change $old @rpath/$new $lib
		fi
		let k=k+1
	done < /tmp/t.lis
done < /tmp/l.lis
# Set links to the new lib
ln -s libXgmt.5.dylib libXgmt.dylib
ln -s libXpsl.5.dylib libXpsl.dylib
# stuff for gs:
cp /opt/local/lib/libgs.9.16.dylib libXgs.9.16.dylib 
cp /opt/local/lib/libfreetype.6.dylib libXfreetype.6.dylib
install_name_tool -id @rpath/libXgs.9.16.dylib libXgs.9.16.dylib 
install_name_tool -id @rpath/libXfreetype.6.dylib libXfreetype.6.dylib
install_name_tool -change /opt/local/lib/libtiff.5.dylib @rpath/libXtiff.5.dylib libXgs.9.16.dylib 
install_name_tool -change /opt/local/lib/libfreetype.6.dylib @rpath/libXfreetype.6.dylib libXgs.9.16.dylib 

# Do plugin supplement separately since not called lib*
echo "--> mex/gmt/plugins/supplements.so"
echo
cd gmt/plugins
otool -L supplements.so | grep executable_path | awk '{print $1}' > /tmp/t.lis
let k=1
while read old; do
	new=`echo $old | awk -F/ '{printf "libX%s\n", substr($NF,4)}'`
	install_name_tool -change $old @rpath/$new supplements.so
	let k=k+1
done < /tmp/t.lis
