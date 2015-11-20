function run_tests(what)
% Run 'packages' of examples/tests
% WHAT can be 'all_examples' to run all examples or 'all_tests' for the ported tests

	if (strcmp(what, 'all_examples'))
		for (k = 1:45)
			ex = sprintf('ex%.2d', k);
			gmtest(ex)
		end
	elseif (strcmp(what, 'all_tests'))
		do_tests()
	end

% -------------------------------------------------------------------
function do_tests()
% Run all tests that we have ported so far
	tests = {
		'poldecimate','gmtspatial'
		'measure','gmtspatial'
		'readwrite_withgdal','grdimage'
		'vectors','grdvector'
		'gspline_5','greenspline'
		'select','ogr'
		'spotter_6','spotter'
		'map_JE','psbasemap'
		'mapscales','psbasemap'
		'oblique','psbasemap'		% THIS ONE CRASHES gmt(['psbasemap -R-100/100/-60/60 -JOa1/0/45/5.5i -B30g30 -P -K -Xc > lixo.ps'])
		'oblsuite','pscoast'
		'trimvector','psxy'
		'sph_5','sph'
		%'flexure_e','potential'	% See test. We dont' have syntax to run it
		'clipping6','psxy'
		%'clipping5','psxy'			% See test. We dont' have syntax to run it
		'clipping4','psxy'
		'geosegmentize','psxy'
		};
	
	for (k = 1:size(tests,1))
		gmtest(tests{k,1}, tests{k,2})
	end