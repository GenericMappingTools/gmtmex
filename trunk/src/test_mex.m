function  test_mex(opt)
%	$Id $
%	Test suite for the GMT-MEX API
%

all_tests = {'gmtread' 'gmtwrite' 'blockmean' 'filter1d' 'psbasemap' 'pscoast' 'gmtinfo' 'surface' 'gmtmath' 'gmtsimplify'}; 

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
			case 'blockmean',	blockmean
			case 'filter1d',	filter1d
			case 'psbasemap',	psbasemap
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

function G = blockmean()
	disp ('Test blockmean');
	gmt('create')
	t = rand(100,3) * 100;
	gmt('blockmean -R0/150/0/100 -I10 > ave1.txt', t);
	B = gmt('blockmean -R0/150/0/100 -I10', t);
	gmt ('write -Td ave2.txt', B);
	gmt('destroy')

function G = filter1d
	disp ('Test filter1d');
	gmt('create')
	t = zeros (100,2);
	t(:,1) = 1:100;
	t(:,2) = rand(100,1) * 100;
	gmt('filter1d -Fg10 -E > filt1.txt', t);
	F = gmt('filter1d -Fg10 -E', t);
	gmt ('write -Td filt2.txt', F);
	gmt('destroy')

function G = surface()
	disp ('Test surface');
	gmt('create')
	t = rand(100,3) * 100;
	% Make a somewhat recognizible grid instead of random numbers
	t(:,3) = (t(:,1)/150).^2 - (t(:,2)/100).^2 + 1.0;
	gmt('surface -R0/150/0/100 -I1 -Glixo.grd', t);
	G = gmt('surface -R0/150/0/100 -I1', t);
	G = gmt('surface -R0/150/0/100 -I1 -G', t);
	gmt('destroy')

function pscoast()
	disp ('Test pscoast');
	gmt('create')
	gmt('pscoast -R110/140/20/35 -JB125/20/25/45/5i -Bag -Dl -Ggreen -Wthinnest -A250 -P > GMT_albers.ps')
	gmt('destroy')

function psbasemap()
	disp ('Test psbasemap');
	gmt('create')
	gmt('psbasemap -R110/140/20/35 -JB125/20/25/45/5i -Bafg -BWSne+ggreen -P > plot.ps')
	gmt('destroy')

function gmtinfo()
	disp ('Test gmtinfo');
	gmt('create')
	t = rand(100,3) * 100;
	r = gmt('info -C', t);
	disp(['gmtinfo of random 0-100 is ' num2str(r)])
	gmt('destroy')

function gmtmath()
	disp ('Test gmtmath');
	gmt('create');
	t1 = rand(10,1);
	t2 = rand(10,1);
	gmt ('write -Td t1.txt', t1);
	gmt ('write -Td t2.txt', t2);
	A = gmt('math t1.txt t2.txt ADD 0.5 MUL LOG10 =');
	if (~isempty(A))
		disp('gmtmath gave something when operating on two input files given on command line')
	end
	B = gmt('math $ $ ADD 0.5 MUL LOG10 =', t1, t2);
	if (~isempty(B))
		disp('gmtmath gave something when operating on two input files given as vectors')
	end
	gmt('destroy')

function gmtsimplify()
	disp ('Test gmtsimplify');
	gmt('create')
	t = rand(50,2);
	t2 = gmt('simplify -T0.2', t);
	gmt('destroy')
