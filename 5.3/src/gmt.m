function varargout = gmt(cmd, varargin)
% Helper function to call the gmtmex MEX function
%

%	$Id: $

	if (nargin == 0)
		disp('Usage: gmt(''module_name options'', numeric_opts)')
		return
	end
	
	[varargout{1:nargout}] = gmtmex(cmd, varargin{:});
