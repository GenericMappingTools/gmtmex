function varargout = gmt(cmd, varargin)
% Helper function to call the gmtmex MEX function

%	$Id$

	if (nargin == 0)
		fprintf(sprintf('Usage: to call a GMT program\n\tgmt(''module_name options'', numeric_opts)\n\n'))
		fprintf(sprintf(['       To create a Grid struct from a 2d Z array and a 1x9 header vector\n\t' ...
		                 'G = gmt(''_fill_grid_struct'', Z, head)\n\n']))
		fprintf(sprintf(['       To create an Image struct from a 2d img array and a 1x9 header vector\n\t' ...
		                 'I = gmt(''_fill_img_struct'', img, head [,cmap])\n' ...
						 '       Here, and above, HEAD is a vector with [x_min x_max, y_min y_max z_min z_max reg x_inc y_inc]\n' ...
						 '       and CMAP is a color palette structure or a Matlab Mx3 cmap array (not yet).\n\n']))
		fprintf(sprintf(['       To join two color palette structures\n\t' ...
		                 'cpt = gmt(''_cptjoin'', cpt1, cpt2)\n']))
		return
	end

	if (cmd(1) == '_')
		[varargout{1:nargout}] = feval(cmd(2:end), varargin{:});
	else
		[varargout{1:nargout}] = gmtmex(cmd, varargin{:});
	end

% -------------------------------------------------------------------------------------------------
function cpt = cptjoin(cpt1, cpt2)
% Join two CPT1 and CPT2 color palette structures. 
% Note, the two palettes should be continuous across its common border. No testing on that is donne here

	if (nargin ~= 2)
		error('    Must provide 2 input arguments.')
	elseif (cpt1.depth ~= cpt2.depth)
		error('    Cannot join two palletes that have different bit depths.')
	end
	if (size(cpt1.colormap,1) ~= size(cpt1.range))
		% A continuous palette so the join would have one color in excess. We could average
		% the top cpt1 color and bottom cpt2 but that would blur the transition. 
		%cpt.colormap = [cpt1.colormap(1:end-1,:); (cpt1.colormap(end,:)+cpt2.colormap(1,:))/2; cpt2.colormap(2:end,:)];
		cpt.colormap = [cpt1.colormap(1:end-1,:); cpt2.colormap];
		cpt.alpha    = [cpt1.alpha(1:end-1,:);    cpt2.alpha];
	else
		cpt.colormap = [cpt1.colormap; cpt2.colormap];
		cpt.alpha    = [cpt1.alpha;    cpt2.alpha];
	end
	cpt.range       = [cpt1.range;    cpt2.range];
	cpt.rangeMinMax = [cpt1.rangeMinMax(1) cpt2.rangeMinMax(2)];
	cpt.BFN         = cpt1.BFN;			% Just keep the first one
	cpt.depth       = cpt1.depth;

% -------------------------------------------------------------------------------------------------
function G = fill_grid_struct(Z, head)
% Fill the Grid struct used in gmtmex. HEAD is the old 1x9 header vector.

	if (nargin ~= 2)
		error('    Must provide 2 input arguments.')
	elseif (size(Z,1) < 2 || size(Z,2) < 2)
		error('    First argin must be a decent 2D array.')
	elseif (any(size(head) ~= [1 9]))
		error('    Second argin must be a 1x9 header vector.')
	end

	if (~isa(head, 'double')),	head = double(head);	end
	G.ProjectionRefPROJ4 = '';
	G.ProjectionRefWKT = '';	
	G.range = head(1:6);
	G.inc = head(8:9);
	G.n_rows = size(Z,1);
	G.n_columns = size(Z,2);
	G.n_bands = size(Z,3);
	G.registration = head(7);
	G.NoDataValue = NaN;
	G.title = '';
	G.remark = '';
	G.command = '';
	G.DataType = 'float32';
	G.x = linspace(head(1), head(2), G.n_columns);
	G.y = linspace(head(3), head(4), G.n_rows);
	G.z = Z;
	G.x_units = '';
	G.y_units = '';
	G.z_units = '';	

% -------------------------------------------------------------------------------------------------
function I = fill_img_struct(img, head, cmap)
% Fill the Image struct used in gmtmex. HEAD is the old 1x9 header vector.

	if (nargin < 2)
		error('    Must provide at least 2 input arguments.')
	end
	if (size(img,1) < 2 || size(img,2) < 2)
		error('    First argin must be a decent 2D image array.')
	elseif (any(size(head) ~= [1 9]))
		error('    Second argin must be a 1x9 header vector.')
	end

	if (~isa(head, 'double')),	head = double(head);	end
	I.ProjectionRefPROJ4 = '';
	I.ProjectionRefWKT = '';	
	I.range = head(1:6);
	I.inc = head(8:9);
	I.n_rows = size(img,1);
	I.n_columns = size(img,2);
	I.n_bands = size(img,3);
	I.NoDataValue = NaN;
	I.registration = head(7);
	I.title = '';
	I.remark = '';
	I.command = '';
	I.DataType = 'uint8';
	I.x = linspace(head(1), head(2), I.n_columns);
	I.y = linspace(head(3), head(4), I.n_rows);
	I.image = img;
	I.x_units = '';
	I.y_units = '';
	I.z_units = '';	
	if (nargin == 3)
		if (~isa(cmap, 'struct'))
			% TODO: write a function that converts from Mx3 Matlab cmap to color struct used in MEX
			error('The third argin (cmap) must be a colormap struct.')
		end
		I.colormap = cmap;	
	else
		I.colormap = [];	
	end
	if (I.n_bands == 4)			% Not obvious that this is the best choice
		I.alpha = img(:,:,4);
		I.n_bands = 3;
	else
		I.alpha = [];	
	end
