function out = gmtest(test, test_dir)
% ...
	global g_root_dir out_path;
	% Edit those two for your own needs
	g_root_dir = 'C:/progs_cygw/GMTdev/gmt5/branches/5.2.0/';
	out_path = 'V:/';		% Set this if you want to save the PS files in a prticular place
	GRAPHICSMAGICK = 'C:/programs/GraphicsMagick/gm.exe';

	if ((nargin == 1) && (numel(test) == 4) && strncmp(test, 'ex', 2))	% Run example from gallery
		[ps, t_path] = gallery(test);
	else
		[ps, t_path] = call_test(test, test_dir);
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
	ps_orig  = [t_path fname '.ps'];

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
function [ps, t_path] = call_test(test, test_dir)
	global g_root_dir out_path
	pato = [g_root_dir 'test/' test_dir];
	addpath(pato)		% Need this because that's the only way I found to make str2func() work
	try
		[ps, t_path] = feval(str2func(test), out_path);
	catch
		disp(lasterr)
		ps = [];	t_path = [];	% Defaults for the case of error
	end
	rmpath(pato)

