# GMT/MEX - GMT API for MATLAB

# Extended Docs

https://github.com/GenericMappingTools/gmtmex/wiki

# Introduction

The GMT MATLAB interface makes it possible to access all GMT modules from MATLAB. Users of MATLAB can write MATLAB scripts
that call upon GMT modules to do any of the things GMT normally can do, and return the results (grids, data-tables, CPTs,
text-files, and even final images via psconvert) to MATLAB variables. MATLAB matrices can be given as input to GMT modules.
Examples below will give you the general idea.

# Installing

## Windows

The Windows installers come already with the gmtmex.mexw64|32 and gmt.m files necessary run the MEX. Only make sure that the
GMT6.4 binary dir is either in the Windows path (the installer does that for you) and in the MATLAB path (you have to do it
yourself). If you want to (re)build the MEX file yourself, see the compile_mex.bat in the source SVN repository.

## macOS

We have successfully built the MATLAB interface under macOS. However, due to the way MATLAB handles shared libraries it is a
delicate process, with several caveats. This may change over time as we work with MathWorks to straighten out the kinks.
The following works:

 * Install the GMT macOS Bundle
 * Run the gmt_prepmex.sh script in the bundle's share/tools directory.  This will duplicate
   the GMT 6.4 installation into /opt/gmt and re-baptize all the shared libraries.
 * Use gmtswitch to make /opt/gmt the current active GMT version
 * Checkout the gmtmex project via git into some directory, i.e.,
   ```
   git clone https://github.com/GenericMappingTools/gmtmex.git
   ```
 * In `gmtmex/`, run autoconf then configure --enable-matlab (and maybe --enable-debug) is you
   can help debug things.
 * Run make which builds the gmtmex.mexmaci64 (x86_64) or gmtmex.mexmaca64 (arm64).  This executable
   is accessed by the gmt.m script.
 * Set your MATLAB path so these two can be found (or copy them to a suitable directory).
 * Make sure your gmt.conf file has the entry GMT_CUSTOM_LIBS=/opt/gmt/lib/gmt/plugins/supplements.so.

You can also build your own bundle (see CMakeLists.txt in main GMT directory).  The above works
with UNIX installations from fink or HomeBrew but fails for us if under MacPorts (then, MATLAB
will complain about wrong shared HDF5 library and we crash).
If you wish to help debug in XCode then see the gmtmex wiki for more details.  While the latest
2023b MATLAB version works with XCode 14, earlier versions may require an older Xcode version.
We used the 2023b MATLAB version to build the interface.

## Unix/Linux

Preliminary experiments indicate we will have to fight the shared library dilemma here as well.
Volunteers on Linux wishing to run the GMT MATLAB interface are needed to make progress.

# Using

The MATLAB wrapper was designed to work in a way the closest as possible to the command line version
and yet to provide all the facilities of the MATLAB IDE (the ML command line desktop). In this sense,
all **GMT** options are put in a single text string that is passed, plus the data itself when it applies,
to the ``gmt()`` command. For example to reproduce the CookBook example of an Hemisphere map using a
Azimuthal projection

    gmt('pscoast -Rg -JA280/30/3.5i -Bg -Dc -A1000 -Gnavy -P > GMT_lambert_az_hemi.ps')

but that is not particularly interesting as after all we could do the exact same thing on the a shell
command line. Things start to get interesting when we can send data *in* and *out* from MATLAB to
**GMT**. So, consider the following example

    t = rand(100,3) * 150;
    G = gmt('surface -R0/150/0/150 -I1', t);

Here we just created a random data *100x3* matrix and told **GMT** to grid it using it's program
*surface*. Note how the syntax follows closely the standard usage but we sent the data to be
interpolated (the *t* matrix) as the second argument to the ``gmt()`` function. And on return we
got the *G* variable that is a structure holding the grid and it's metadata. See the
``grid struct`` for the details of its members.

Imagining that we want to plot that random data art, we can do it with a call to *grdimage*, like

    gmt('grdimage -JX8c -Ba -P -Cblue,red > crap_img.ps', G)

Note that we now sent the *G grid* as argument instead of the **-G** *gridname* that we would have
used in the command line. But for readability we could well had left the **-G** option in command string. E.g:

    gmt('grdimage -JX8c -Ba -P -Cblue,red -G > crap_img.ps', G)

While for this particular case it makes no difference to use or not the **-G**, because there is **only**
one input, the same does not hold true when we have more than one. For example, we can run the same example
but compute the CPT separately.

    cpt = gmt('grd2cpt -Cblue,red', G);
    gmt('grdimage -JX8c -Ba -P -C -G > crap_img.ps', G, cpt)

Now we had to explicitly write the **-C** & **-G** (well, actually we could have omitted the **-G** because
it's a mandatory input but that would make the things more confusing). Note also the order of the input data variables.
It is crucial that any *required* (primary) input data objects (for grdimage that is the grid) are given before
any *optional* (secondary) input data objects (here, that is the CPT object).  The same is true for modules that
return more than one item: List the required output object first followed by optional ones.

To illustrate another aspect on the importance of the order of input data let us see how to plot a sine curve
made of colored filled circles.

    x = linspace(-pi, pi)';            % The *xx* var
    seno = sin(x);                     % *yy*
    xyz  = [x seno seno];              % Duplicate *yy* so that it can be colored
    cpt  = gmt('makecpt -T-1/1/0.1');  % Create a CPT
    gmt('psxy -R-3.2/3.2/-1.1/1.1 -JX12c -Sc0.1c -C -P -Ba > seno.ps', xyz, cpt)

The point here is that we had to give *xyz, cpt* and not *cpt, xyz* (which would error) because optional input data
associated with an option letter **always comes after the required input**.

To plot text strings we send in the input data wrapped in a cell array. Example:

    lines = {'5 6 Some label', '6 7 Another label'};
    gmt('pstext -R0/10/0/10 -JM6i -Bafg -F+f18p -P > text.ps', lines)

and we get back text info in cell arrays as well. Using the *G* grid computed above we can run *gmtinfo* on it

    info = gmt('info', G)

At the end of an **GMT** session work we call the internal functions that will do the house keeping of
freeing no longer needed memory. We do that with this command:

    gmt('destroy')


So that's basically how it works. When numeric data have to be sent *in* to **GMT** we use
MATLAB variables holding the data in matrices or structures or cell arrays, depending on data type. On
return we get the computed result stored in variables that we gave as output arguments.
Things only complicate a little more for the cases where we can have more than one *input* or
*output* arguments, since the order or the arguments matter (Remember the rule: primary first, secondary second).
The file *gallery.m*, that reproduces the examples in the Gallery section of the GMT
documentation, has many (not so trivial) examples on usage of the MEX GMT API.
