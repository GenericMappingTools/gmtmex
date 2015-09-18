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
		'vectors','gmtvector'
		};
	
	for (k = 1:size(tests,1))
		gmtest(tests{k,1}, tests{k,2})
	end