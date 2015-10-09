#!/bin/bash
#
# Until/if we are able to get Matlab to not override every
# single request for a shared library with its own out-of-date
# version which results in version conflicts, we have to use
# this trick under OS X:
#
# Duplicate the lib, bin, include files in /opt/gmt
# Rebaptize libs with unique names by inserting "X"
# Link the gmt.mex executable with these libraries.
#
# To prepare your system to run the gmt.mex application, run
# /Application/GMT-5.2.x[_r#####.app]/Contents/Resources/share/tools/gmt_prepmex.sh
#
# First get a reliable absolute path to the bundle's top directory
pushd `dirname $0` > /dev/null
BUNDLEDIR=`pwd | sed -e sB/Contents/Resources/share/toolsBBg`
popd > /dev/null
# Set path to our new gmt installation
MEXGM5TDIR=/opt/gmt
# Set path to additional subdirectories
MEXLIBDIR=$MEXGM5TDIR/lib
MEXINCDIR=$MEXGM5TDIR/include
MEXBINDIR=$MEXGM5TDIR/bin
MEXSUPDIR=$MEXLIBDIR/gmt/plugins
# Create install directory [remove first if exist]
sudo rm -rf $MEXGM5TDIR
sudo mkdir -p $MEXBINDIR $MEXSUPDIR $MEXINCDIR
# Find user's group and use that to set ownership
grp=`id -gn`
sudo chown -R ${USER}:${grp} $MEXGM5TDIR
# Place the include files
cd $BUNDLEDIR/Contents/Resources/include
scp -r gmt $MEXINCDIR
# Place the bin files
cd $BUNDLEDIR/Contents/Resources/bin
scp -r * $MEXBINDIR
# Now do lib files
cd $BUNDLEDIR/Contents/Resources/lib
# Find a list of all libs shipped with the OSX bundle, except ours:
ls *.dylib | egrep -v 'libgmt.dylib|libpostscriptlight.dylib' > /tmp/l.lis
# For each, duplicate into /opt/gmt but add a leading X to each name
while read lib; do
	new=`echo $lib | awk '{printf "libX%s\n", substr($1,4)}'`
	cp $lib $MEXLIBDIR/$new
done < /tmp/l.lis
# Same for the supplement shared plugin
cp gmt/plugins/supplements.so $MEXLIBDIR/gmt/plugins/Xsupplements.so
cd $MEXLIBDIR
ls *.dylib > /tmp/l.lis
# For all libs in $MEXLIBDIR, change internal references to contain the leading "X"
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
			install_name_tool -id $MEXLIBDIR/$new $lib
		else
			install_name_tool -change $old $MEXLIBDIR/$new $lib
		fi
		let k=k+1
	done < /tmp/t.lis
done < /tmp/l.lis
# Set links to the new libs
ln -s libXgmt.dylib libgmt.dylib
ln -s libXgmt.5.dylib libXgmt.dylib
ln -s libXpostscriptlight.5.dylib libXpostscriptlight.dylib
# Same stuff for gs which is called by psconvert as a system call:
cp /opt/local/lib/libgs.9.16.dylib libXgs.9.16.dylib 
cp /opt/local/lib/libfreetype.6.dylib libXfreetype.6.dylib
install_name_tool -id $MEXLIBDIR/libXgs.9.16.dylib libXgs.9.16.dylib 
install_name_tool -id $MEXLIBDIR/libXfreetype.6.dylib libXfreetype.6.dylib
install_name_tool -change /opt/local/lib/libtiff.5.dylib $MEXLIBDIR/libXtiff.5.dylib libXgs.9.16.dylib 
install_name_tool -change /opt/local/lib/libfreetype.6.dylib $MEXLIBDIR/libXfreetype.6.dylib libXgs.9.16.dylib 

# Do plugin supplement separately since not called lib*
echo "--> gmt/plugins/supplements.so"
echo
cd gmt/plugins
otool -L Xsupplements.so | grep executable_path | awk '{print $1}' > /tmp/t.lis
let k=1
while read old; do
	new=`echo $old | awk -F/ '{printf "libX%s\n", substr($NF,4)}'`
	install_name_tool -change $old $MEXLIBDIR/$new Xsupplements.so
	let k=k+1
done < /tmp/t.lis

# Do bin dir
cd $MEXBINDIR
otool -L gmt | grep executable_path | awk '{print $1}' > /tmp/t.lis
let k=1
while read old; do
	new=`echo $old | awk -F/ '{printf "libX%s\n", substr($NF,4)}'`
	install_name_tool -change $old $MEXLIBDIR/$new gmt
	let k=k+1
done < /tmp/t.lis

# Finally fix gmt-config
cat << EOF > /tmp/skip
GMT_EXEDIR=
CONFIG_CFLAGS=
CONFIG_INCLUDEDIR=
CONFIG_LIBS=
CONFIG_PREFIX=
EOF
sed '/GMT_EXEDIR/q' gmt-config > /tmp/new
cat << EOF >> /tmp/new
CONFIG_CFLAGS="-I/opt/gmt/include/gmt"
CONFIG_DATA=\$(\$GMT_EXEDIR/gmt --show-datadir)
CONFIG_INCLUDEDIR="/opt/gmt/include/gmt"
CONFIG_LIBS="-L/opt/gmt/lib -lXgmt"
CONFIG_PREFIX="/opt/gmt"
EOF
sed -n '/GMT_EXEDIR/,$p' gmt-config | grep -v -f/tmp/skip >> /tmp/new
mv -f /tmp/new gmt-config
chmod +x gmt-config
