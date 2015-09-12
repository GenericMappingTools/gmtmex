function run_tests(what)
%	$Id$
% Run 'packages' of examples/tests
	if (strcmp(what, 'all_examples'))
		for (k = 1:45)
			ex = sprintf('ex%.2d', k);
			gmtest(ex)
		end
	else
	end

