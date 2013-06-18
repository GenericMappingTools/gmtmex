function  test_mex(opt)

all_tests = {'gmtread' 'minmax' 'pscoast' 'surface'}; 

if (nargin == 0)
	opt = all_tests;
else
	opt = {opt};		% Make it a cell to fit the other branch
end


try
	for (k = 1: numel(opt))
		switch opt{k}
			case 'gmtread'		% 
				gmtread;
			case 'minmax'
				minmax;
			case 'pscoast'
				pscoast
			case 'surface'
				surface;
		end
	end
catch
	sprintf('Error in test: %s\n%s', opt{k}, lasterr)
end

function G = gmtread()
	surface;			% Create the lixo.grd grid
	gmt('create')
	G = gmt('gmtread -Tg lixo.grd');
	gmt('grdcontour -JX6i -P -Ba -C5 -Vd ->crap.ps', G);
	gmt('destroy')

function G = surface()
	t = rand(100,3) * 100;
	gmt('create')
	gmt('surface -R0/150/0/100 -I1 -Glixo.grd', t);
	G = gmt('surface -R0/150/0/100 -I1', t);
	%G = gmt('surface -R0/150/0/100 -I1 -G', t);		% This crashes
	gmt('destroy')

function pscoast()
	gmt('create')
	gmt('pscoast -R110/140/20/35 -JB125/20/25/45/5i -Bag -Dl -Ggreen -Wthinnest -A250 -P > GMT_albers.ps')
	gmt('destroy')

function minmax()
	t = rand(100,3) * 100;
	gmt('create')
	r = gmt('minmax $', t);
	r = gmt('minmax', t);
	disp(['minmax of random 0-100 is ' r])
	gmt('destroy')

