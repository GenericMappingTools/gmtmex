function  gallery(opt)
%	$Id$
%	The examples Gallery in GMT-MEX API
%

global g_root_dir;	g_root_dir = 'C:/progs_cygw/GMTdev/gmt5/branches/5.2.0/';
global out_path;	out_path = 'V:/';		% Set this if you want to save the PS files in a prticular place

	all_exs = {'ex01' 'ex02' 'ex23' 'ex44'}; 

	if (nargin == 0)
		opt = all_exs;
	else
		opt = {opt};		% Make it a cell to fit the other branch
	end

	try
		for (k = 1: numel(opt))
			switch opt{k}
				case 'ex01',   ex01
				case 'ex02',   ex02
				case 'ex23',   ex23
				case 'ex44',   ex44
			end
		end
	catch
		sprintf('Error in test: %s\n%s', opt{k}, lasterr)
	end

	gmt('destroy')

% -------------------------------------------------------------------------------------------------
function ex01()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex01'];
	ps = [out_path 'example_01.ps'];

	gmt('gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0 FONT_ANNOT_PRIMARY 10p')
	gmt(['psbasemap -R0/6.5/0/7.5 -Jx1i -B0 -P -K > ' ps])
	gmt(['pscoast -Rg -JH0/6i -X0.25i -Y0.2i -O -K -Bg30 -Dc -Glightbrown -Slightblue >> ' ps])
	cmd = sprintf('grdcontour %s/osu91a1f_16.nc', d_path);
	gmt([cmd ' -J -C10 -A50+f7p -Gd4i -L-1000/-1 -Wcthinnest,- -Wathin,- -O -K -T+d0.1i/0.02i >> ' ps])
	gmt([cmd ' -J -C10 -A50+f7p -Gd4i -L-1/1000 -O -K -T+d0.1i/0.02i >> ' ps])
	gmt(['pscoast -Rg -JH6i -Y3.4i -O -K -B+t"Low Order Geoid" -Bg30 -Dc -Glightbrown -Slightblue >> ' ps])
	gmt([cmd ' -J -C10 -A50+f7p -Gd4i -L-1000/-1 -Wcthinnest,- -Wathin,- -O -K -T+d0.1i/0.02i+l >> ' ps])
	gmt([cmd ' -J -C10 -A50+f7p -Gd4i -L-1/1000 -O -T+d0.1i/0.02i+l >> ' ps])
	builtin('delete','gmt.conf');

% -------------------------------------------------------------------------------------------------
function ex02()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex02'];
	ps = [out_path 'example_02.ps'];

	gmt('gmtset FONT_TITLE 30p MAP_ANNOT_OBLIQUE 0')
	g_cpt = gmt('makecpt -Crainbow -T-2/14/2');
	cmd = sprintf('grdimage %s/HI_geoid2.nc', d_path);
	gmt([cmd ' -R160/20/220/30r -JOc190/25.5/292/69/4.5i -E50 -K -P -B10 -X1.5i -Y1.25i > '  ps], g_cpt)
	gmt(['psscale -DjRM+o0.6i/0+jLM+w2.88i/0.4i+mc+e -R -J -O -K -Bx2+lGEOID -By+lm >> ' ps], g_cpt)
	t_cpt = gmt(sprintf('grd2cpt %s/HI_topo2.nc -Crelief -Z', d_path));
	GHI_topo2_int = gmt(sprintf('grdgradient %s/HI_topo2.nc -A0 -Nt', d_path));
	cmd = sprintf('grdimage %s/HI_topo2.nc', d_path);
	gmt([cmd ' -I -R -J -B+t"H@#awaiian@# T@#opo and @#G@#eoid@#"' ...
        ' -B10 -E50 -O -K -C -Y4.5i --MAP_TITLE_OFFSET=0.5i >> ' ps], GHI_topo2_int, t_cpt)
	gmt(['psscale -DjRM+o0.6i/0+jLM+w2.88i/0.4i+mc -R -J -O -K -I0.3 -Bx2+lTOPO -By+lkm >> ' ps], t_cpt)
% gmt pstext -R0/8.5/0/11 -Jx1i -F+f30p,Helvetica-Bold+jCB -O -N -Y-4.5i >> $ps << END
% -0.4 7.5 a)
% -0.4 3.0 b)
% END
	builtin('delete','gmt.conf');

% -------------------------------------------------------------------------------------------------
function ex23()
	global out_path
	ps = [out_path 'example_23.ps'];

	% Position and name of central point:
	lon  = 12.50;	lat  = 41.99;
	name = 'Rome';

	Gdist = gmt(sprintf('grdmath -Rg -I1 %f %f SDIST', lon, lat));

	gmt(['pscoast -Rg -JH90/9i -Glightgreen -Sblue -A1000 -Dc -Bg30 -B+t"Distances from ' ...
		name ' to the World" -K -Wthinnest > ' ps])

	gmt(['grdcontour -A1000+v+u" km"+fwhite -Glz-/z+ -S8 -C500 -O -K -JH90/9i -Wathin,white ' ...
		' -Wcthinnest,white,- >> ' ps], Gdist)
	
	% Location info for 5 other cities + label justification
	cities = cell(5,1);
	cities{1} = '105.87 21.02 LM HANOI';
	cities{2} = '282.95 -12.1 LM LIMA';
	cities{3} = '178.42 -18.13 LM SUVA';
	cities{4} = '237.67 47.58 RM SEATTLE';
	cities{5} = '28.20 -25.75 LM PRETORIA';
	fid = fopen('cities.d','w');
	for (k = 1:5)
		fprintf(fid, '%s\n', cities{k});
	end
	fclose(fid);

	% For each of the cities, plot great circle arc to Rome with gmt psxy
	gmt(['psxy -R -J -O -K -Wthickest,red >> ' ps], [lon lat; 105.87 21.02])
	gmt(['psxy -R -J -O -K -Wthickest,red >> ' ps], [lon lat; 282.95 -12.1])
	gmt(['psxy -R -J -O -K -Wthickest,red >> ' ps], [lon lat; 178.42 -18.13])
	gmt(['psxy -R -J -O -K -Wthickest,red >> ' ps], [lon lat; 237.67 47.58])
	gmt(['psxy -R -J -O -K -Wthickest,red >> ' ps], [lon lat; 28.20 -25.75])

	% Plot red squares at cities and plot names:
	for (k = 1:5)
		gmt(['pstext -R -J -O -K -Dj0.15/0 -F+f12p,Courier-Bold,red+j -N >> ' ps], cities(k))	% It should accept a plain string as well
	end

	% Place a yellow star at Rome
	gmt(['psxy -R -J -O -K -Sa0.2i -Gyellow -Wthin >> ' ps], [12.5 41.99])

	% Sample the distance grid at the cities and use the distance in km for labels
	dist = gmt('grdtrack cities.d -G', Gdist);
	t = cell(5,1);
	for (k = 1:5)
		t{k} = sprintf('%f %f %d', dist(k,1), dist(k,2), round(dist(k,end)));
	end
	gmt(['pstext -R -J -O -D0/-0.2i -N -Gwhite -W -C0.02i -F+f12p,Helvetica-Bold+jCT >> ' ps], t)
	builtin('delete','cities.d');

% -------------------------------------------------------------------------------------------------
function ex44()
	global out_path
	ps = [out_path 'example_44.ps'];

	% Bottom map of Australia
	gmt(['pscoast -R110E/170E/44S/9S -JM6i -P -Baf -BWSne -Wfaint -N2/1p  -EAU+gbisque -Gbrown' ...
		' -Sazure1 -Da -K -Xc --FORMAT_GEO_MAP=dddF > ' ps])
	gmt(['psbasemap -R -J -O -K -DjTR+w1.5i+o0.15i/0.1i+sxx000 -F+gwhite+p1p+c0.1c+s >> ' ps])
	t = load('xx000');		% x0 y0 w h
	cmd = sprintf('pscoast -Rg -JG120/30S/%f -Da -Gbrown -A5000 -Bg -Wfaint -EAU+gbisque -O -K -X%f -Y%f >> %s', t(3), t(1), t(2), ps);
	gmt(cmd)
	gmt(sprintf('psxy -R -J -O -K -T -X-%f -Y-%f >> %s', t(1), t(2), ps))
	% Determine size of insert map of Europe
	t = gmt('mapproject -R15W/35E/30N/48N -JM2i -W');	% w h
	gmt(['pscoast -R10W/5E/35N/44N -JM6i -Baf -BWSne -EES+gbisque -Gbrown -Wfaint -N1/1p -Sazure1' ...
		' -Df -O -K -Y4.5i --FORMAT_GEO_MAP=dddF >> ' ps])
	gmt(sprintf('psbasemap -R -J -O -K -DjTR+w%f/%f+o0.15i/0.1i+sxx000 -F+gwhite+p1p+c0.1c+s >> %s', t(1), t(2), ps))
	t = load('xx000');		% x0 y0 w h
	cmd = sprintf('pscoast -R15W/35E/30N/48N -JM%f -Da -Gbrown -B0 -EES+gbisque -O -K -X%f -Y%f ', t(3), t(1), t(2));
	gmt([cmd '--MAP_FRAME_TYPE=plain >> ' ps])
	gmt(sprintf('psxy -R -J -O -T -X-%f -Y-%f >> %s', t(1), t(2), ps))
	builtin('delete','xx000');
