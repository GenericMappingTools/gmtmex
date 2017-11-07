function run_tests(what)
% Run 'packages' of examples/tests
% WHAT can be 'all_examples' to run all examples, 'all_tests' for the ported tests
%      or a single example. e.g. run_tests('ex02')

	if (strcmp(what, 'all_examples'))
		for (k = 1:46)
			ex = sprintf('ex%.2d', k);
			gmtest(ex)
		end
	elseif (strncmp(what, 'ex', 2))
		gmtest(what)
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
		'oblique','psbasemap'
		'oblsuite','pscoast'
		'trimvector','psxy'
		'decoratedlines','psxy'
		'sph_5','sph'
		%'flexure_e','potential'	% See test. We dont' have syntax to run it
		'spheres','potential'
		'clipping6','psxy'
		%'clipping5','psxy'			% See test. We dont' have syntax to run it
		'clipping4','psxy'
		'geosegmentize','psxy'
		};
	
	for (k = 1:size(tests,1))
		gmtest(tests{k,1}, tests{k,2})
	end