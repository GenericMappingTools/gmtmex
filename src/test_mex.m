function  test_mex(opt)
%	$Id$
%	Test suite for the GMT-MEX API
%

all_tests = {'gmtread' 'gmtinfo' 'pscoast' 'surface' 'gmtmath' 'simplify'}; 

if (nargin == 0)
	opt = all_tests;
else
	opt = {opt};		% Make it a cell to fit the other branch
end

try
	for (k = 1: numel(opt))
		switch opt{k}
			case 'gmtread',		gmtread;
			case 'gmtwrite',	gmtwrite;
			case 'pscoast',		pscoast
			case 'gmtinfo',		gmtinfo;
			case 'surface',		surface;
			case 'gmtmath',		gmtmath;
			case 'gmtsimplify',	gmtsimplify;
		end
	end
catch
	sprintf('Error in test: %s\n%s', opt{k}, lasterr)
	gmt('destroy')
end

function G = gmtread()
	disp ('Test read');
	surface;			% Create the lixo.grd grid
	gmt('create')
	G = gmt('read -Tg lixo.grd');
	gmt('grdcontour -JX6i -P -Ba -C5 > crap.ps', G);
	gmt('destroy')

function G = gmtwrite()
	disp ('Test write');
	gmt('create')
	G = gmt('read -Tg lixo.grd');
	gmt ('write -Tg crap.grd', G);
	gmt('destroy')

function G = surface()
	disp ('Test surface');
	t = rand(100,3) * 100;
	gmt('create')
	gmt('surface -R0/150/0/100 -I1 -Glixo.grd', t);
	G = gmt('surface -R0/150/0/100 -I1', t);
	G = gmt('surface -R0/150/0/100 -I1 -G', t);
	gmt('destroy')

function pscoast()
	disp ('Test pscoast');
	gmt('create')
	gmt('pscoast -R110/140/20/35 -JB125/20/25/45/5i -Bag -Dl -Ggreen -Wthinnest -A250 -P > GMT_albers.ps')
	gmt('destroy')

function minmax()
	disp ('Test gmtinfo');
	t = rand(100,3) * 100;
	gmt('create')
	r = gmt('info -C', t);
	disp(['gmtinfo of random 0-100 is ' num2str(r)])
	gmt('destroy')

function gmtmath()
	disp ('Test gmtmath');
	t1 = rand(10,1);
	t2 = rand(10,1);
	gmt('create')
	t = gmt('math ADD 0.5 MUL LOG10 =', t1, t2);
	if (~isempty(t))
		disp('gmtmath gave something when operating into two input files')
	end
	gmt('destroy')

function gmtsimplify()
	disp ('Test gmtsimplify');
	t = rand(50,2);
	gmt('create')
	t2 = gmt('simplify -T0.01', t);
	gmt('destroy')
