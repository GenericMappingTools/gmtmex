function  test_mex(opt)
%	$Id$
%	Test suite for the GMT-MEX API
%

all_tests = {'blockmean' 'filter1d' 'gmtinfo' 'gmtmath' 'gmtread' 'gmtsimplify' 'gmtwrite' 'psbasemap' ...
	'pscoast' 'pstext' 'psxy' 'grd2xyz' 'grdinfo' 'grdimage' 'grdsample' 'surface', 'coasts'}; 

if (nargin == 0)
	opt = all_tests;
else
	opt = {opt};		% Make it a cell to fit the other branch
end

try
	for (k = 1: numel(opt))
		switch opt{k}
			case 'blockmean',   blockmean
			case 'filter1d',    filter1d
			case 'gmtinfo',    	gmtinfo;
			case 'gmtmath',     gmtmath;
			case 'gmtread',     gmtread;
			case 'gmtsimplify',	gmtsimplify;
			case 'gmtwrite',    gmtwrite;
			case 'psbasemap',   psbasemap
			case 'pscoast',    	pscoast
			case 'pstext',    	pstext
			case 'psxy',    	psxy
			case 'grd2xyz',     grd2xyz;
			case 'grdinfo',    	grdinfo;
			case 'grdimage',    grdimage;
			case 'grdsample',   grdsample;
			case 'surface',     surface;
			case 'coasts',      coasts;
		end
	end
catch
	sprintf('Error in test: %s\n%s', opt{k}, lasterr)
	gmt('destroy')
end

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
	C = gmt('math $ $ ADD 0.5 MUL LOG10', t1, t2);
	if (~isempty(C))
		disp('gmtmath gave something when operating on two input files given as vectors')
	end
	gmt('destroy')

function G = gmtread()
	disp ('Test read');
	surface;			% Create the lixo.grd grid
	gmt('create')
	G = gmt('read -Tg lixo.grd');
	gmt('grdcontour -JX6i -P -Ba -C0.2 > crap.ps', G);
	gmt('destroy')

function gmtsimplify()
	disp ('Test gmtsimplify');
	gmt('create')
	t = rand(50,2);
	t2 = gmt('simplify -T0.2', t);
	gmt('destroy')

function G = gmtwrite()
	disp ('Test write');
	gmt('create')
	G = gmt('read -Tg lixo.grd');
	gmt ('write -Tg crap.grd', G);
	gmt('destroy')

function grd2xyz()
	disp ('Test grd2xyz');
	gmt('create')
	G = gmt('surface -R0/150/0/150 -I1', rand(100,3) * 100);
	xyz = gmt('grd2xyz', G);
	gmt('destroy')

function grdinfo()
	disp ('Test grdinfo');
	gmt('create')
	G = gmt('surface -R0/150/0/150 -I1', rand(100,3) * 100);
	T = gmt('grdinfo', G);
	gmt('destroy')

function grdimage()
	disp ('Test grdimage');
	gmt('create')
	t = rand(100,3) * 100;
	G = gmt('surface -R0/150/0/100 -I1', t);
	gmt('grdimage -JX8c -Ba -P -Cblue,red > crap_img.ps', G);
	gmt('destroy')

function grdsample()
	disp ('Test grdsample');
	gmt('create')
	G = gmt('surface -R0/150/0/150 -I1', rand(100,3) * 100);
	gmt('grdsample -I100+/100+', G);
	gmt('destroy')

function psbasemap()
	disp ('Test psbasemap');
	gmt('create')
	gmt('psbasemap -R110/140/20/35 -JB125/20/25/45/5i -Bafg -BWSne+ggreen -P > plot.ps')
	gmt('destroy')

function pscoast()
	disp ('Test pscoast');
	gmt('create')
	gmt('pscoast -R110/140/20/35 -JB125/20/25/45/5i -Bag -Dl -Ggreen -Wthinnest -A250 -P > GMT_albers.ps')
	gmt('destroy')

function pstext()
	disp ('Test pstext');
	gmt('create')
	lines = {'5 6 Some label', '6 7 Another label'};
	gmt('pstext -R0/10/0/10 -JM6i -Bafg -F+f18p -P > text.ps', lines)
	gmt('destroy')

function psxy()
	disp ('Test psxy');
	gmt('create')
	gmt('psxy -Sc0.5c -G191/101/95  -JX10c -R0/10/0/10 -W1,171/43/33  -P -K >  Jlogo.ps', [0.5 1.0])
	gmt('psxy -Sc0.4c -G158/122/190 -JX10c -R0/10/0/10 -W1,130/83/171 -O -K >> Jlogo.ps', [1.5 1.0])
	gmt('psxy -Sc0.4c -G128/171/93  -JX10c -R0/10/0/10 -W1,81/143/24  -O    >> Jlogo.ps', [1.0 1.5])
	gmt('destroy')

function grdcontour()
	disp ('Test grdcontour after surface');
	gmt('create')
	G = gmt('read -Tg lixo.grd');
	D = gmt('grdcontour -C0.1 -D', G);
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

function coasts()
	disp ('Test coastlines');
	opt_res = ' -Di';	opt_N = ' -Na';		opt_I = ' -Ia';	opt_R = ' -R-10/10/30/50';
	gmt('create')
	coast = gmt(['pscoast -M -W ' opt_R, opt_res, ' -A1/1/1']);
	boundaries = gmt(['pscoast -M ' opt_R, opt_N, opt_res]);
	rivers = gmt(['pscoast -M ', opt_R, opt_I, opt_res]);
	gmt('destroy')
	h = figure; hold on
	plot(coast(:,1), coast(:,2))
	plot(boundaries(:,1), boundaries(:,2))
	plot(rivers(:,1), rivers(:,2))
	hold off
	%pause(1.0);		delete(h);
