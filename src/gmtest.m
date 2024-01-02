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

	global GMT_ROOT_DIR GMT_PLOT_DIR GMT_GM_EXE
	% Edit those two for your own needs
	if (~exist('GMT_ROOT_DIR', 'var') || isempty(GMT_ROOT_DIR))
		GMT_ROOT_DIR = 'C:/progs_cygw/GMTdev/gmt5/trunk/';
		GMT_PLOT_DIR = 'V:/';		% Set this if you want to save the PS files in a prticular place
		GMT_GM_EXE  = 'C:/programs/GraphicsMagick/gm.exe';
	end

	if (nargin < 3),	family = '';	end

	if ((nargin == 1) && (numel(test) == 4) && strncmp(test, 'ex', 2))	% Run example from gallery
		[ps, orig_path] = gmt_gallery(test, GMT_ROOT_DIR, GMT_PLOT_DIR);
	else
		[ps, orig_path] = call_test(test, test_dir, GMT_ROOT_DIR, GMT_PLOT_DIR, family);
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
	cmd = [GMT_GM_EXE ' compare -density 200 -maximum-error 0.005 -highlight-color magenta -highlight-style' ...
		' assign -metric rmse -file ' png_name ' ' ps_orig ' ' ps];
	[status, cmdout] = system(cmd);

	n_lines = find(cmdout == sprintf('\n'));
 	if (status == 1)		% GM returned an error status
		if (numel(n_lines) == 8)
			t = cmdout(n_lines(end-1)+1:end-1);		% Last line containing the failing info
			ind = strfind(t, 'image');
			disp(['Test ' test ' FAIL with: ' t(ind:end)])
		else
			disp(['GMT_GM_EXE failed to run test: ' test])
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
			orig_path = orig_path(1:end-3);		% Because original PS live in a subdir below
		end
	catch
		disp(lasterr)
	end
	rmpath(pato)

