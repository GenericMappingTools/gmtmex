function out = gmtest(test)
	GRAPHICSMAGICK = 'C:/programs/GraphicsMagick/gm.exe';

	[ps, t_path] = gallery(test);
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
	
% 	if (status == 1)
		n_lines = find(cmdout == sprintf('\n'));
		if (numel(n_lines) == 8)
			t = cmdout(n_lines(end-1)+1:end-1);		% Last line containing the failing info
			ind = strfind(t, 'image');
			disp(['Test ' test ' FAIL with: ' t(ind:end)])
		else
			disp(['Test ' test ' PASS'])
			builtin('delete', png_name);
		end
% 	else
% 		disp(['GRAPHICSMAGICK failed to run test: ' test])
% 	end
	
	if (nargout)
		out = cmdout;
	end
	
	