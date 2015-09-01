function  gallery(opt)
%	$Id$
%	The examples Gallery in GMT-MEX API
%

global g_root_dir out_path;
% Edit those two for your own needs
g_root_dir = 'C:/progs_cygw/GMTdev/gmt5/branches/5.2.0/';
out_path = 'V:/';		% Set this if you want to save the PS files in a prticular place

	all_exs = {'ex01' 'ex02' 'ex04' 'ex06' 'ex07' 'ex08' 'ex09' 'ex10' 'ex12' 'ex13' 'ex14' 'ex15' ...
		'ex23' 'ex34' 'ex35' 'ex36' 'ex38' 'ex39' 'ex40' 'ex41' 'ex44'}; 

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
				case 'ex04',   ex04
				case 'ex06',   ex06
				case 'ex07',   ex07
				case 'ex08',   ex08
				case 'ex09',   ex09
				case 'ex10',   ex10
				case 'ex12',   ex12
				case 'ex13',   ex13
				case 'ex14',   ex14
				case 'ex15',   ex15
				case 'ex23',   ex23
				case 'ex34',   ex34
				case 'ex35',   ex35
				case 'ex36',   ex36
				case 'ex38',   ex38
				case 'ex39',   ex39
				case 'ex40',   ex40
				case 'ex41',   ex41
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
function ex04()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex04'];
	ps = [out_path 'example_04.ps'];

	fid = fopen('zero.cpt','w');
	fprintf(fid, '%s\n', '-10  255   0  255');
	fprintf(fid, '%s\n', '  0  100  10  100');
	fclose(fid);

	cmd = sprintf('grdcontour %s/HI_geoid4.nc', d_path);
	gmt([cmd ' -R195/210/18/25 -Jm0.45i -p60/30 -C1 -A5+o -Gd4i -K -P -X1.25i -Y1.25i > ' ps])
	gmt(['pscoast -R -J -p -B2 -BNEsw -Gblack -O -K -TdjBR+o0.1i+w1i+l >> ' ps])
	gmt([sprintf('grdview %s/HI_topo4.nc', d_path) ' -R195/210/18/25/-6/4 -J -Jz0.34i -p -Czero.cpt -O -K ' ...
		' -N-6+glightgray -Qsm -B2 -Bz2+l"Topo (km)" -BneswZ -Y2.2i >> ' ps])
	gmt(['pstext -R0/10/0/10 -Jx1i -F+f60p,ZapfChancery-MediumItalic+jCB -O >> ' ps], {'3.25 5.75 H@#awaiian@# R@#idge@#'})
	builtin('delete','zero.cpt');

	ps = [out_path 'example_04c.ps'];
	Gg_intens = gmt([sprintf('grdgradient %s/HI_geoid4.nc', d_path) ' -A0 -Nt0.75 -fg']);
	Gt_intens = gmt([sprintf('grdgradient %s/HI_topo4.nc', d_path)  ' -A0 -Nt0.75 -fg']);
	gmt([sprintf('grdimage %s/HI_geoid4.nc', d_path) ...
		' -I -R195/210/18/25 -JM6.75i -p60/30 -C' d_path '/geoid.cpt -E100 -K -P -X1.25i -Y1.25i > ' ps], Gg_intens)
	gmt(['pscoast -R -J -p -B2 -BNEsw -Gblack -O -K >> ' ps])
	gmt(['psbasemap -R -J -p -O -K -TdjBR+o0.1i+w1i+l --COLOR_BACKGROUND=red --FONT=red' ...
		' --MAP_TICK_PEN_PRIMARY=thinner,red >> ' ps])
	gmt(['psscale -R -J -p240/30 -DJBC+o0/0.5i+w5i/0.3i+h -C' d_path '/geoid.cpt -I -O -K -Bx2+l"Geoid (m)" >> ' ps])
	gmt([sprintf('grdview %s/HI_topo4.nc', d_path) ' -I -R195/210/18/25/-6/4 -J -C' d_path '/topo.cpt' ...
		' -JZ3.4i -p60/30 -O -K -N-6+glightgray -Qc100 -B2 -Bz2+l"Topo (km)" -BneswZ -Y2.2i >> ' ps], Gt_intens)
	gmt(['pstext -R0/10/0/10 -Jx1i -F+f60p,ZapfChancery-MediumItalic+jCB -O >> ' ps], {'3.25 5.75 H@#awaiian@# R@#idge@#'})

% -------------------------------------------------------------------------------------------------
function ex06()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex06'];
	ps = [out_path 'example_06.ps'];

	gmt(['psrose ' d_path '/fractures.d -: -A10r -S1.8in -P -Gorange -R0/1/0/360 -X2.5i -K -Bx0.2g0.2' ...
		' -By30g30 -B+glightblue -W1p > ' ps])
	gmt(['pshistogram -Bxa2000f1000+l"Topography (m)" -Bya10f5+l"Frequency"+u" %"' ...
		' -BWSne+t"Histograms"+glightblue ' d_path '/v3206.t -R-6000/0/0/30 -JX4.8i/2.4i -Gorange -O' ...
		' -Y5.0i -X-0.5i -L1p -Z1 -W250 >> ' ps])

% -------------------------------------------------------------------------------------------------
function ex07()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex07'];
	ps = [out_path 'example_07.ps'];

	gmt(['pscoast -R-50/0/-10/20 -JM9i -K -Slightblue -GP300/26:FtanBdarkbrown -Dl -Wthinnest' ...
		' -B10 --FORMAT_GEO_MAP=dddF > ' ps])
	gmt(['psxy -R -J -O -K ' d_path '/fz.xy -Wthinner,- >> ' ps])
	gmt(['psxy ' d_path '/quakes.xym -R -J -O -K -h1 -Sci -i0,1,2s0.01 -Gred -Wthinnest >> ' ps])
	gmt(['psxy -R -J -O -K ' d_path '/isochron.xy -Wthin,blue >> ' ps])
	gmt(['psxy -R -J -O -K ' d_path '/ridge.xy -Wthicker,orange >> ' ps])
	gmt(['psxy -R -J -O -K -Gwhite -Wthick -A >> ' ps], [-14.5 15.2; -2 15.2; -2 17.8; -14.5 17.8])
	gmt(['psxy -R -J -O -K -Gwhite -Wthinner -A >> ' ps], [-14.35 15.35; -2.15 15.35; -2.15 17.65; -14.35 17.65])
	gmt(['psxy -R -J -O -K -Sc0.08i -Gred -Wthinner >> ' ps], [-13.5 16.5])
	gmt(['pstext -R -J -F+f18p,Times-Italic+jLM -O -K >> ' ps], {'-12.5 16.5 ISC Earthquakes'})
	gmt(['pstext -R -J -O -F+f30,Helvetica-Bold,white=thin >> ' ps], {'-43 -5 SOUTH' '-43 -8 AMERICA' '-7 11 AFRICA'})

% -------------------------------------------------------------------------------------------------
function ex08()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex08'];
	ps = [out_path 'example_08.ps'];

	xyz = gmt(['grd2xyz ' d_path '/guinea_bay.nc']);
	gmt(['psxyz -B1 -Bz1000+l"Topography (m)" -BWSneZ+b+tETOPO5' ...
		' -R-0.1/5.1/-0.1/5.1/-5000/0 -JM5i -JZ6i -p200/30 -So0.0833333ub-5000 -P' ...
		' -Wthinnest -Glightgreen -K > ' ps], xyz)
	gmt(['pstext -R -J -JZ -Z0 -F+f24p,Helvetica-Bold+jTL -p -O >> ' ps], {'0.1 4.9 This is the surface of cube'})

% -------------------------------------------------------------------------------------------------
function ex09()
% THIS EXAMPLE FAILS
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex09'];
	ps = [out_path 'example_09.ps'];

	gmt(['pswiggle ' d_path '/tracks.txt -R185/250/-68/-42 -K -Jm0.13i -Ba10f5 -BWSne+g240/255/240 -G+red' ...
		' -G-blue -Z2000 -Wthinnest -S240/-67/500/@~m@~rad --FORMAT_GEO_MAP=dddF > ' ps])
	gmt(['psxy -R -J -O -K ' d_path '/ridge.xy -Wthicker >> ' ps])
	gmt(['psxy -R -J -O -K ' d_path '/fz.xy    -Wthinner,- >> ' ps])
	% Take label from segment header and plot near coordinates of last record of each track
	% BUT THIS CURRENTLY FAILS BECAUSE CODE FLOW GOES TRHOUGH gmt_api #2492 AND 
	% GMT_IS_REFERENCE_VIA_VECTOR DOESN'T KNOW HOW TO DEAL WITH THE gmtconvert -El OPTION
	%t = gmt(['gmtconvert -El ' d_path '/tracks.txt']);
	%gmt(['pstext -R -J -F+f10p,Helvetica-Bold+a50+jRM+h -D-0.05i/-0.05i -O >> ' ps], t)

% -------------------------------------------------------------------------------------------------
function ex10()
% THIS EXAMPLE FAILS
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex10'];
	ps = [out_path 'example_10.ps'];

	gmt(['pscoast -Rd -JX8id/5id -Dc -Sazure2 -Gwheat -Wfaint -A5000 -p200/40 -K > ' ps])
	fid = fopen([d_path '/languages.txt']);
	str = fread(fid,'*char');
	fclose(fid);
	str = strread(str','%s','delimiter','\n');
	k = 1;
	while (str{k}(1) == '#')
		k = k + 1;
	end
	str = str(k:end);		% Remove the comment lines
	nl = numel(str);
	array = zeros(nl, 7);
	for (k = 1:nl)
		array(k,:) = strread(str{k}, '%f', 7);
	end
	t = cell(nl,1);
	for (k = 1:nl)
		t{k} = sprintf('%d %d %d\n',array(k,1:2),sum(array(k,3:end)));
	end
	gmt(['pstext -R -J -O -K -p -Gwhite@30 -D-0.25i/0 -F+f30p,Helvetica-Bold,firebrick=thinner+jRM >> ' ps], t)
	gmt(['psxyz ' d_path '/languages.txt -R-180/180/-90/90/0/2500 -J -JZ2.5i -So0.3i -Gpurple -Wthinner' ...
		' --FONT_TITLE=30p,Times-Bold --MAP_TITLE_OFFSET=-0.7i -O -K -p --FORMAT_GEO_MAP=dddF' ...
		' -Bx60 -By30 -Bza500+lLanguages -BWSneZ+t"World Languages By Continent" >> ' ps])
	gmt(['psxyz -R -J -JZ -So0.3ib -Gblue -Wthinner -O -K -p >> ' ps], [array(:,1:2) sum(array(:,3:4),2) array(:,3)])
	gmt(['psxyz -R -J -JZ -So0.3ib -Gdarkgreen -Wthinner -O -K -p >> ' ps], [array(:,1:2) sum(array(:,3:5),2) sum(array(:,3:4),2)])
	gmt(['psxyz -R -J -JZ -So0.3ib -Gyellow -Wthinner -O -K -p >> ' ps], [array(:,1:2) sum(array(:,3:6),2) sum(array(:,3:5),2)])
	gmt(['psxyz -R -J -JZ -So0.3ib -Gred -Wthinner -O -K -p >> ' ps], [array(:,1:2) sum(array(:,3:7),2) sum(array(:,3:6),2)])
	% BUG HERE. NOTHING IS WRITTEN BY THE NEXT COMMAND 
	gmt(['pslegend -R -J -JZ -DjLB+o0.2i+w1.35i/0+jBL -O --FONT=Helvetica-Bold' ...
		' -F+glightgrey+pthinner+s-4p/-6p/grey20@40 -p ' d_path '/legend.txt -Vl >> ' ps])

% -------------------------------------------------------------------------------------------------
function ex12()
% THIS EXAMPLE FAILS
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex12'];
	ps = [out_path 'example_12.ps'];

	net_xy = gmt(['triangulate ' d_path '/table_5.11 -M']);
	gmt(['psxy -R0/6.5/-0.2/6.5 -JX3.06i/3.15i -B2f1 -BWSNe -Wthinner -P -K -X0.9i -Y4.65i > ' ps], net_xy)
	gmt(['psxy ' d_path '/table_5.11 -R -J -O -K -Sc0.12i -Gwhite -Wthinnest >> ' ps])
	t = load([d_path '/table_5.11']);
	nl = size(t,1);
	c = cell(nl,1);
	for (k = 1:nl)
		c{k} = sprintf('%f %f %d\n', t(k,1:2), k-1);
	end
	gmt(['pstext -R -J -F+f6p -O -K >> ' ps], c)
	%
	% Then draw network and print the node values
	%
	gmt(['psxy -R -J -B2f1 -BeSNw -Wthinner -O -K -X3.25i >> ' ps], net_xy)
	gmt(['psxy -R -J -O -K ' d_path '/table_5.11 -Sc0.03i -Gblack >> ' ps])
	gmt(['pstext ' d_path '/table_5.11 -R -J -F+f6p+jLM -O -K -Gwhite -W -C0.01i -D0.08i/0i -N >> ' ps])
	%
	% Then contour the data and draw triangles using dashed pen; use "gmt gmtinfo" and "gmt makecpt" to make a
	% color palette (.cpt) file
	%
	T = gmt(['info -T25/2 ' d_path '/table_5.11']);
	topo_cpt = gmt(['makecpt -Cjet ' T{1}]);
	gmt(['pscontour -R -J ' d_path '/table_5.11 -B2f1 -BWSne -Wthin -C -Lthinnest,-' ...
		' -Gd1i -X-3.25i -Y-3.65i -O -K >> ' ps], topo_cpt)
	%
	% Finally color the topography
	% NEX CALL ERRORS with "pscontour: -I option requires constant color between contours!"		BECAUSE
% in GMT_Read_Data #L4786
% (void)GMTAPI_split_via_method (API, API->object[item]->method, &via); 
% sets via = 0, which ends in
% GMT_init_cpt #L2998
% P->is_continuous = true;	/* The only kind from the outside (?) */
	gmt(['pscontour -R -J ' d_path '/table_5.11 -B2f1 -BeSnw -C -I -X3.25i -O -K >> ' ps], topo_cpt)
	gmt(['pstext -R0/8/0/11 -Jx1i -F+f30p,Helvetica-Bold+jCB -O -X-3.25i >> ' ps], {'3.16 8 Delaunay Triangulation'})

% -------------------------------------------------------------------------------------------------
function ex13()
% THIS EXAMPLE FAILS
	global out_path
	ps = [out_path 'example_13.ps'];

	Gz = gmt('grdmath -R-2/2/-2/2 -I0.1 X Y R2 NEG EXP X MUL =');
	Gdzdx = gmt('grdmath $ DDX', Gz);
	Gdzdy = gmt('grdmath $ DDY', Gz);
	gmt(['grdcontour -JX3i -B1 -BWSne -C0.1 -A0.5 -K -P -Gd2i -S4 -T+d0.1i/0.03i > ' ps], Gdzdx)
	gmt(['grdcontour -J -B -C0.05 -A0.2 -O -K -Gd2i -S4 -T+d0.1i/0.03i -Xa3.45i >> ' ps], Gdzdy)
	gmt(['grdcontour -J -B -C0.05 -A0.1 -O -K -Gd2i -S4 -T+d0.1i/0.03i -Y3.45i >> ' ps], Gz)
	gmt(['grdcontour -J -B -C0.05 -O -K -Gd2i -S4 -X3.45i >> ' ps], Gz)
	gmt(['grdvector  $ $ -I0.2 -J -O -K -Q0.1i+e+n0.25i -Gblack -W1p -S5i --MAP_VECTOR_SHAPE=0.5 >> ' ps], Gdzdx, Gdzdy)
	gmt(['pstext -R0/6/0/4.5 -Jx1i -F+f40p,Times-Italic+jCB -O -X-3.45i >> ' ps], ...
		{'3.2 3.6 z(x,y) = x@~\327@~exp(-x@+2@+-y@+2@+)'})

% -------------------------------------------------------------------------------------------------
function ex14()
% THIS EXAMPLE FAILS
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex14'];
	ps = [out_path 'example_14.ps'];

	% First draw network and label the nodes
	gmt('gmtset MAP_GRID_PEN_PRIMARY thinnest,-')
	gmt(['psxy ' d_path '/table_5.11 -R0/7/0/7 -JX3.06i/3.15i -B2f1 -BWSNe -Sc0.05i -Gblack -P -K -Y6.45i > ' ps])
	gmt(['pstext ' d_path '/table_5.11 -R -J -D0.1c/0 -F+f6p+jLM -O -K -N >> ' ps])
	mean_xyz = gmt(['blockmean ' d_path '/table_5.11 -R0/7/0/7 -I1']);

	% Then draw gmt blockmean cells
	gmt(['psbasemap -R0.5/7.5/0.5/7.5 -J -O -K -Bg1 -X3.25i >> ' ps])
	gmt(['psxy -R0/7/0/7 -J -B2f1 -BeSNw -Ss0.05i -Gblack -O -K >> ' ps], mean_xyz)
	% Reformat to one decimal for annotation purposes
%	t = gmt('gmtconvert --FORMAT_FLOAT_OUT=%.1f', mean_xyz);							% <---------------- FAILS HERE
%	gmt(['pstext -R -J -D0.15c/0 -F+f6p+jLM -O -K -Gwhite -W -C0.01i -N >> ' ps], t)

	% Then gmt surface and contour the data
	Gdata = gmt('surface -R -I1', mean_xyz);
	gmt(['grdcontour -J -B2f1 -BWSne -C25 -A50 -Gd3i -S4 -O -K -X-3.25i -Y-3.55i >> ' ps], Gdata)
	gmt(['psxy -R -J -Ss0.05i -Gblack -O -K >> ' ps], mean_xyz)

	% Fit bicubic trend to data and compare to gridded gmt surface
	Gtrend = gmt('grdtrend -N10 -T', Gdata);
	track  = gmt('project -C0/0 -E7/7 -G0.1 -N');
	gmt(['grdcontour -J -B2f1 -BwSne -C25 -A50 -Glct/cb -S4 -O -K -X3.25i >> ' ps], Gtrend)
	gmt(['psxy -R -J -Wthick,. -O -K >> ' ps], track)

	% Sample along diagonal
	data  = gmt('grdtrack -G -o2,3', Gdata, track);
	trend = gmt('grdtrack -G -o2,3', Gtrend, track);
	t = gmt('info -I0.5/25', data, trend);
	gmt(['psxy -JX6.3i/1.4i ' t{1} ' -Wthick -O -K -X-3.25i -Y-1.9i -Bx1 -By50 -BWSne >> ' ps], data)
	gmt(['psxy -R -J -Wthinner,- -O >> ' ps], trend)

% -------------------------------------------------------------------------------------------------
function ex15()
% THIS EXAMPLE FAILS TO PLOT THE STAR AT THE MINIMUM AT UR FIG
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex15'];
	ps = [out_path 'example_15.ps'];

	gmt(['gmtconvert ' d_path '/ship.xyz -bo > ship.b'])

	region = gmt('info ship.b -I1 -bi3d');
	region = region{1};				% We want this to be astring, not a cell
	Gship  = gmt(['nearneighbor ' region ' -I10m -S40k ship.b -bi']);
	gmt(['grdcontour -JM3i -P -B2 -BWSne -C250 -A1000 -Gd2i -K > ' ps], Gship)

	ship_10m = gmt(['blockmedian ' region ' -I10m ship.b -b3d']);
	Gship = gmt(['surface ' region ' -I10m'], ship_10m);
	gmt(['psmask ' region ' -I10m ship.b -J -O -K -T -Glightgray -bi3d -X3.6i >> ' ps])
	gmt(['grdcontour -J -B -C250 -L-8000/0 -A1000 -Gd2i -O -K >> ' ps], Gship)

	gmt(['psmask ' region ' -I10m -J -B -O -K -X-3.6i -Y3.75i >> ' ps], ship_10m)
	gmt(['grdcontour -J -C250 -A1000 -L-8000/0 -Gd2i -O -K >> ' ps], Gship)
	gmt(['psmask -C -O -K >> ' ps])

	Gship_clipped = gmt('grdclip -Sa-1/NaN -G', Gship);
	gmt(['grdcontour -J -B -C250 -A1000 -L-8000/0 -Gd2i -O -K -X3.6i >> ' ps], Gship_clipped)
	gmt(['pscoast ' region ' -J -O -K -Ggray -Wthinnest >> ' ps])
	info = gmt('grdinfo -C -M', Gship);
	info = strread(info{1}, '%f', 14);
	gmt(['psxy -R -J -O -K -Sa0.15i -Wthick >> ' ps], info(11:12))		% <--------- DOES NOT SHOW UP
	gmt(['pstext -R0/3/0/4 -Jx1i -F+f24p,Helvetica-Bold+jCB -O -N >> ' ps], {'-0.3 3.6 Gridding with missing data'})
	builtin('delete','ship.b');

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
function ex34()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex34'];
	ps = [out_path 'example_34.ps'];

	gmt('gmtset FORMAT_GEO_MAP dddF')
	gmt(['pscoast -JM4.5i -R-6/20/35/52 -EFR,IT+gP300/8 -Glightgray -Baf -BWSne -P -K -X2i > ' ps])
	% Extract a subset of ETOPO2m for this part of Europe
	% gmt grdcut etopo2m_grd.nc -R -GFR+IT.nc=ns
	z_cpt = gmt('makecpt -Cglobe -T-5000/5000/500 -Z');
	FR_IT_int = gmt(['grdgradient ' d_path '/FR+IT.nc -A15 -Ne0.75 -G']);
	gmt(['grdimage ' d_path '/FR+IT.nc -I -C -J -O -K -Y4.5i' ...
		' -Baf -BWsnE+t"Franco-Italian Union, 2042-45" >> ' ps], FR_IT_int, z_cpt)	% Hmmm, how does it know which input is which?
	gmt(['pscoast -J -R -EFR,IT+gred@60 -O >> ' ps])
	builtin('delete','gmt.conf');

% -------------------------------------------------------------------------------------------------
function ex35()
% THIS EXAMPLE FAILS BECAUSE OF sphtriangulate
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex35'];
	ps = [out_path 'example_35.ps'];

	% Get the crude GSHHS data, select GMT format, and decimate to ~20%:
	% gshhs $GMTHOME/src/coast/gshhs/gshhs_c.b | $AWK '{if ($1 == ">" || NR%5 == 0) print $0}' > gshhs_c.txt
	% Get Voronoi polygons
	tt_pol = gmt(['sphtriangulate ' d_path '/gshhs_c.txt -Qv -D']);
	% Compute distances in km
	Gtt = gmt('sphdistance -Rg -I1 -Q -G -Lk', tt_pol);
	t_cpt = gmt('makecpt -Chot -T0/3500/500 -Z');
	% Make a basic image plot and overlay contours, Voronoi polygons and coastlines
	gmt(['grdimage -JG-140/30/7i -P -K -C -X0.75i -Y2i > ' ps], Gtt, t_cpt)
	gmt(['grdcontour -J -O -K -C500 -A1000+f10p,Helvetica,white -L500' ...
		' -GL0/90/203/-10,175/60/170/-30,-50/30/220/-5 -Wa0.75p,white -Wc0.25p,white >> ' ps], Gtt)
	gmt(['psxy -R -J -O -K -W0.25p,green,. >> ' ps], tt_pol)
	gmt(['pscoast -R -J -O -W1p -Gsteelblue -A0/1/1 -B30g30 -B+t"Distances from GSHHG crude coastlines" >> ' ps])
	builtin('delete','gmt.conf');

% -------------------------------------------------------------------------------------------------
function ex36()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex36'];
	ps = [out_path 'example_36.ps'];

	% Interpolate data of Mars radius from Mariner9 and Viking Orbiter spacecrafts
	tt_cpt = gmt('makecpt -Crainbow -T-7000/15000/1000 -Z');
	% Piecewise linear interpolation; no tension
	Gtt = gmt(['sphinterpolate ' d_path '/mars370.txt -Rg -I1 -Q0 -G']);
	gmt(['grdimage -JH0/6i -Bag -C -P -Xc -Y7.25i -K > ' ps], tt_cpt, Gtt)
	gmt(['psxy -Rg -J -O -K ' d_path '/mars370.txt -Sc0.05i -G0 -B30g30 -Y-3.25i >> ' ps])
	% Smoothing
	Gtt = gmt(['sphinterpolate ' d_path '/mars370.txt -Rg -I1 -Q3 -G']);
	gmt(['grdimage -J -Bag -C -Y-3.25i -O -K >> ' ps], tt_cpt, Gtt)
	gmt(['psxy -Rg -J -O -T >> ' ps])

% -------------------------------------------------------------------------------------------------
function ex38()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex38'];
	ps = [out_path 'example_38.ps'];

	t_cpt = gmt('makecpt -Crainbow -T0/1700/100 -Z');
	c_cpt = gmt('makecpt -Crainbow -T0/15/1');
	Gitopo = gmt(['grdgradient ' d_path '/topo.nc -Nt1 -fg -A45 -G']);
	Gout  = gmt(['grdhisteq ' d_path '/topo.nc -G -C16']);
	gmt(['grdimage ' d_path '/topo.nc -I -C -JM3i -Y5i -K -P -B5 -BWSne > ' ps], Gitopo, t_cpt)
	gmt(['pstext -R' d_path '/topo.nc -J -O -K -F+jTR+f14p -T -Gwhite -W1p -Dj0.1i >> ' ps], {'315 -10 Original'})
	gmt(['grdimage -C -J -X3.5i -K -O -B5 -BWSne >> ' ps], c_cpt, Gout)
	gmt(['pstext -R -J -O -K -F+jTR+f14p -T -Gwhite -W1p -Dj0.1i >> ' ps], {'315 -10 Equalized'})
	gmt(['psscale -Dx0i/-0.4i+jTC+w5i/0.15i+h+e+n -O -K -C -Ba500 -By+lm >> ' ps], t_cpt)
	Gout = gmt(['grdhisteq ' d_path '/topo.nc -G -N']);
	c_cpt = gmt('makecpt -Crainbow -T-3/3/0.1 -Z');
	gmt(['grdimage -C -J -X-3.5i -Y-3.3i -K -O -B5 -BWSne >> ' ps], c_cpt, Gout)
	gmt(['pstext -R -J -O -K -F+jTR+f14p -T -Gwhite -W1p -Dj0.1i >> ' ps], {'315 -10 Normalized'})
	Gout = gmt(['grdhisteq ' d_path '/topo.nc -G -N']);
	gmt(['grdimage -C -J -X3.5i -K -O -B5 -BWSne >> ' ps], c_cpt, Gout)
	gmt(['pstext -R -J -O -K -F+jTR+f14p -T -Gwhite -W1p -Dj0.1i >> ' ps], {'315 -10 Quadratic'})
	gmt(['psscale -Dx0i/-0.4i+w5i/0.15i+h+jTC+e+n -O -C -Bx1 -By+lz >> ' ps], c_cpt)

% -------------------------------------------------------------------------------------------------
function ex39()
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex39'];
	ps = [out_path 'example_39.ps'];

	% Evaluate the first 180, 90, and 30 order/degrees of Venus spherical
	% harmonics topography model, skipping the L = 0 term (radial mean).
	% File truncated from http://www.ipgp.fr/~wieczor/SH/VenusTopo180.txt.zip
	% Wieczorek, M. A., Gravity and topography of the terrestrial planets,
	%   Treatise on Geophysics, 10, 165-205, doi:10.1016/B978-044452748-6/00156-5, 2007

	Gv1 = gmt(['sph2grd ' d_path '/VenusTopo180.txt -I1 -Rg -Ng -G -F1/1/25/30']);
	Gv2 = gmt(['sph2grd ' d_path '/VenusTopo180.txt -I1 -Rg -Ng -G -F1/1/85/90']);
	Gv3 = gmt(['sph2grd ' d_path '/VenusTopo180.txt -I1 -Rg -Ng -G -F1/1/170/180']);
	t_cpt = gmt('grd2cpt -Crainbow -E16 -Z', Gv3);
	Gvint = gmt('grdgradient -Nt0.75 -A45 -G', Gv1);
	gmt(['grdimage -I -JG90/30/5i -P -K -Bg -C -X3i -Y1.1i > ' ps], Gvint, t_cpt, Gv1)
	gmt(['pstext -R0/6/0/6 -Jx1i -O -K -Dj0.2i -F+f16p+jLM -N >> ' ps], {'4 4.5 L = 30'})
	gmt(['psscale --FORMAT_FLOAT_MAP="%''g" -C -O -K -Dx1.25i/-0.2i+jTC+w5.5i/0.1i+h -Bxaf -By+lm >> ' ps], t_cpt)
	Gvint = gmt('grdgradient -Nt0.75 -A45 -G', Gv2);
	gmt(['grdimage -I -JG -O -K -Bg -C -X-1.25i -Y1.9i >> ' ps], Gvint, t_cpt, Gv2)
	gmt(['pstext -R0/6/0/6 -Jx1i -O -K -Dj0.2i -F+f16p+jLM -N >> ' ps], {'4 4.5 L = 90'})
	Gv3 = gmt(['sph2grd ' d_path '/VenusTopo180.txt -I1 -Rg -Ng -G -F1/1/170/180']);
	Gvint = gmt('grdgradient -Nt0.75 -A45 -G', Gv3);
	gmt(['grdimage -I -JG -O -K -Bg -C -X-1.25i -Y1.9i >> ' ps], Gvint, t_cpt, Gv3)
	gmt(['pstext -R0/6/0/6 -Jx1i -O -K -Dj0.2i -F+f16p+jLM -N >> ' ps], {'4 4.5 L = 180'})
	gmt(['pstext -R0/6/0/6 -Jx1i -O -F+f24p+jCM -N >> ' ps], {'3.75 5.4 Venus Spherical Harmonic Model'})

% -------------------------------------------------------------------------------------------------
function ex40()
% THIS EXAMPLE FAILS BECAUSE CALLS to gmtspatial CRASH ML
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex40'];
	ps = [out_path 'example_40.ps'];

	% This call crashes ML because:
	%centroid = gmt(['gmtspatial ' d_path '/GSHHS_h_Australia.txt -fg -Qk'])
% GMT_Destroy_Options is called twice. One by gmtspatial and the other at the end of gmt.c (the MEX one).
% In this second call
%  		GMT_free (API->GMT, delete);		/* Then free the structure which was allocated by GMT_memory */ 
% Kabooms because the *head had junk in it. 

	centroid = [133.913549887	-22.9337944115	7592694.55567];
	gmt(['psbasemap -R112/154/-40/-10 -JM5.5i -P -K -B20 -BWSne+g240/255/240 -Xc > ' ps])
	gmt(['psxy ' d_path '/GSHHS_h_Australia.txt -R -J -O -Wfaint -G240/240/255 -K >> ' ps])
	gmt(['psxy ' d_path '/GSHHS_h_Australia.txt -R -J -O -Sc0.01c -Gred -K >> ' ps])
	T500k = gmt(['gmtsimplify ' d_path '/GSHHS_h_Australia.txt -T500k']);
	t = gmt(['gmtspatial ' d_path '/GSHHS_h_Australia.txt -fg -Qk']);
	area = sprintf('Full area = %.0f km@+2@+\n', t(3));
	%| awk '{printf "Full area = %.0f km@+2@+\n", $3}' > area.txt
	t = gmt('gmtspatial -fg -Qk', T500k); 
	area_T500k = sprintf('Reduced area = %.0f km@+2@+\n', t(3));
	%| awk '{printf "Reduced area = %.0f km@+2@+\n", $3}' > area_T500k.txt
	gmt(['psxy -R -J -O -K -W1p,blue >> ' ps], T500k)
	gmt(['psxy -R -J -O -K -Sx0.3i -W3p >> ' ps], centroid)
	gmt(['pstext -R -J -O -K -Dj0.1i/0.1i -F+jTL+f18p >> ' ps], {'112 -10 T = 500 km'})
	gmt(['pstext -R -J -O -K -F+14p+cCM >> ' ps], {area})
	gmt(['pstext -R -J -O -K -F+14p+cLB -Dj0.2i >> ' ps], {area_T500k})
	gmt(['psbasemap -R -J -O -K -B20+lightgray -BWsne+g240/255/240 -Y4.7i >> ' ps])
	gmt(['psxy ' d_path '/GSHHS_h_Australia.txt -R -J -O -Wfaint -G240/240/255 -K >> ' ps])
	gmt(['psxy ' d_path '/GSHHS_h_Australia.txt -R -J -O -Sc0.01c -Gred -K >> ' ps])
	T100k = gmt(['gmtsimplify ' d_path '/GSHHS_h_Australia.txt -T100k']);
	t = gmt('gmtspatial -fg -Qk', T100k');
	area_T100k = sprintf('Reduced area = %.0f km@+2@+\n', t(3));
	%| awk '{printf "Reduced area = %.0f km@+2@+\n", $3}' > area_T100k.txt
	gmt(['psxy -R -J -O -K -W1p,blue >> ' ps], T100k)
	gmt(['psxy -R -J -O -K -Sx0.3i -W3p >> ' ps], centroid)
	gmt(['pstext -R -J -O -K -Dj0.1i/0.1i -F+jTL+f18p >> ' ps], {'112 -10 T = 100 km'})
	gmt(['pstext -R -J -O -K -F+14p+cCM >> ' ps], {area})
	gmt(['pstext -R -J -O -K -F+14p+cLB -Dj0.2i >> ' ps], {area_T100k})
	gmt(['psxy -R -J -O -T >> ' ps])

% -------------------------------------------------------------------------------------------------
function ex41()
% THIS EXAMPLE FAILS, LEGEND IS NOT PLOTTED
	global g_root_dir out_path
	d_path = [g_root_dir 'doc/examples/ex41'];
	ps = [out_path 'example_41.ps'];

	gmt('gmtset FONT_ANNOT_PRIMARY 12p FONT_LABEL 12p')
	gmt(['pscoast -R130W/50W/8N/56N -JM5.6i -B0 -P -K -Glightgray -Sazure1 -A1000 -Wfaint -Xc -Y1.2i --MAP_FRAME_TYPE=plain > ' ps])
	gmt(['pscoast -R -J -O -K -EUS+glightyellow+pfaint -ECU+glightred+pfaint -EMX+glightgreen+pfaint -ECA+glightblue+pfaint >> ' ps])
	gmt(['pscoast -R -J -O -K -N1/1p,darkred -A1000/2/2 -Wfaint -Cazure1 >> ' ps])
	gmt(['psxy -R -J -O -K -Sk' d_path '/my_symbol/0.1i -C' d_path '/my_color.cpt -W0.25p -: ' d_path '/my_data.txt >> ' ps])
	gmt(['pslegend -R0/6/0/9.1 -Jx1i -Dx3i/4.5i+w5.6i+jBC+l1.2 -C0.05i -F+p+gsnow1 -B0 -O ' d_path '/my_table.txt -X-0.2i -Y-0.2i >> ' ps])
	builtin('delete','gmt.conf');

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
