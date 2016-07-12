function out = gmtest(test, test_dir, family)
% Run a  test from the GALLERY, tests or doc/scripts
%
% gmtest(test) run an example from the gallery were TEST is is either 'ex01' ... 'ex45'
%
% gmtest(test, test_dir) runs a test from the tests suit. TEST_DIR is the directory name where TEST lieve
%		example: gmtest('poldecimate','gmtspatial')
%
% gmtest(test, test_dir, family) runs a test from a certainly family of tests.
%	Currently only family = 'scripts' is implemented (means run test from the 'scripts' dir. TEST_DIR is ignored
%		example: gmtest('GMT_insert', '', 'scripts')

%	$Id$

	global g_root_dir out_path;
	% Edit those two for your own needs
	%g_root_dir = '/Users/pwessel/GMTdev/gmt5-dev/trunk/';
	%out_path = '/Users/pwessel/GMTdev/gmt-mex/trunk/src/test/';		% Set this if you want to save the PS files in a prticular place
	%GRAPHICSMAGICK = '/opt/local/bin/gm';
	g_root_dir = 'C:/progs_cygw/GMTdev/gmt5/branches/5.2.0/';
	out_path = 'V:/';		% Set this if you want to save the PS files in a prticular place
	GRAPHICSMAGICK = 'C:/programs/GraphicsMagick/gm.exe';

	if (nargin < 3),	family = '';	end

	if ((nargin == 1) && (numel(test) == 4) && strncmp(test, 'ex', 2))	% Run example from gallery
		[ps, orig_path] = gallery(test, g_root_dir, out_path);
	else
		[ps, orig_path] = call_test(test, test_dir, g_root_dir, out_path, family);
		if (isa(ps, 'logical'))		% In this cases we have no PS files but only an error status from function
			if (ps),		disp(['Test ' test ' PASS'])
			else			disp(['Test ' test ' FAIL'])
			end
			if (nargout),	out = [];	end
			return
		end
	end

	if (isempty(ps))
		disp(['Test -> ' test ' <- does not exist or some other error occurred']),	return
	end
	[pato, fname] = fileparts(ps);
	png_name = [pato filesep fname '.png'];
	ps_orig  = [orig_path fname '.ps'];

	% Compare the ps file with its original.
	cmd = [GRAPHICSMAGICK ' compare -density 200 -maximum-error 0.001 -highlight-color magenta -highlight-style' ...
		' assign -metric rmse -file ' png_name ' ' ps_orig ' ' ps];
	[status, cmdout] = system(cmd);
	
	n_lines = find(cmdout == sprintf('\n'));
 	if (status == 1)		% GM returned an error status
		if (numel(n_lines) == 8)
			t = cmdout(n_lines(end-1)+1:end-1);		% Last line containing the failing info
			ind = strfind(t, 'image');
			disp(['Test ' test ' FAIL with: ' t(ind:end)])
		else
			disp(['GRAPHICSMAGICK failed to run test: ' test])
			disp(cmdout)
		end
	else
		disp(['Test ' test ' PASS'])
		builtin('delete', png_name);
	end
	
	if (nargout)
		out = cmdout;
	end
	
% ---------------------------------------------------------------------------------------
function [ps, orig_path] = call_test(test, test_dir, g_root_dir, out_path, family)
% Here PS hold the full name of the created PS file and ORIG_PATH must
% contain the path to where the original postscript file, against which the PS will be compared

	ps = '';	orig_path = '';		% Defaults for the case of error
	if (isempty(family))						% A test from the tests suit
		pato = [g_root_dir 'test/' test_dir];
	elseif (strcmp(family, 'scripts'))			% A test for the 'scripts' directory
		pato = [g_root_dir 'doc/scripts/ml/'];
	else
		disp(['Family ' family ' is unknown or not yet implemented.'])
		return
	end

	addpath(pato)		% Need this because that's the only way I found to make str2func() work
	try
		[ps, orig_path] = feval(str2func(test), out_path);
		if (strcmp(family, 'scripts'))
			orig_path = orig_path(1:end-3);		% Because original PS lieve in a subdir below
		end
	catch
		disp(lasterr)
	end
	rmpath(pato)

