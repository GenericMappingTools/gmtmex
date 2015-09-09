function out = gmtest(test)
%	$Id$
%	Minimal mimic of the gmtest bash script

	GRAPHICSMAGICK = 'C:/programs/GraphicsMagick/gm.exe';

	[ps, t_path] = gallery(test);
	[pato, fname] = fileparts(ps);
	png_name = [pato filesep fname '.png'];
	ps_orig  = [t_path fname '.ps'];
% 	cmd = ['gmt(psconvert -Tf -A -P ' psname ')'];
% 	system(cmd)

	% Compare the ps file with its original. Check $ps against original $ps or against $1.ps (if $1 given)
	
	% syntax: gm compare [ options ... ] reference-image [ options ... ] compare-image [ options ... ]
	%rms=$("${GRAPHICSMAGICK}" compare -density 200 -maximum-error 0.001 -highlight-color magenta -highlight-style assign -metric rmse -file "${ps%.ps}.png" "$ps" "$src/${psref:-$ps}") || pscmpfailed="yes"
	cmd = [GRAPHICSMAGICK ' compare -density 200 -maximum-error 0.001 -highlight-color magenta -highlight-style' ...
		' assign -metric rmse -file ' png_name ' ' ps_orig ' ' ps];
	[status, cmdout] = system(cmd);
	
	if (status == 1)
		n_lines = find(cmdout == sprintf('\n'));
		if (numel(n_lines == 8))
			t = cmdout(n_lines(end-1)+1:end-1);		% Last line containing the failing info
			ind = strfind(t, 'image');
			disp(['Test ' test ' FAIL with: ' t(ind:end)])
		else
			disp(['Test ' test ' PASS'])
		end
	else
		disp(['GRAPHICSMAGICK failed to run test: ' test])
	end
	
	if (nargout)
		out = cmdout;
	end
	
	
