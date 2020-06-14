function [matrixCoerced, varargout] = Coerce(matrix, varargin)
  % [matrix, changed] = Coerce(matrix, minimum, maximum)
  % Coerce matrix into range defined by minimum and maximum. By default
  % minimum = 0, maximum = 1.
  % changed is true if any value had to be coerced
  %
  %
  % This script and its functions follow the coding style that can be
  % sumarized in:
  % * Variables have lower camel case
  % * Functions upper camel case
  % * Constants all upper case
  % * Spaces around operators
  %
  % Authors:  Néstor Uribe-Patarroyo
  %
  % NUP:
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA;
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  % MGH Therapy monitoring project (v1.0)
  %
  % Changelog:
  %
  % V1.0 (2017-11-01): Initial version released
  %
  % Copyright Néstor Uribe-Patarroyo (2017)

  minimum = 0;
  maximum = 1;
  if nargin >= 2
    minimum = varargin{1};
  end
  if nargin >= 3
    maximum = varargin{2};
  end
  matrixCoerced = max(matrix, minimum);
  matrixCoerced = min(matrixCoerced, maximum);
  
  if nargout > 1
    changed = any(matrix(:) ~= matrixCoerced(:));
    varargout{1} = changed;
  end
end