function [VI, T] = scatteredInterpn(X, V, XI, varargin)
%% SCATTEREDINTERPN Linear interpolation of N-dimensional scattered data
%
% Syntax:
%  VI = scatteredInterpn(X, V, XI)
%  VI = scatteredInterpn(X, V, XI, T)
%  [VI, T] = scatteredInterpn(X, V, XI)
%
%  VI = scatteredInterpn(___, method)
%  VI = scatteredInterpn(___, method, extrapval)
%
% Input:
%  X   - coordinates of points         - [nPoints, nDim]
%  V   - value(s) at points            - [nPoints, nData]
%  XI  - interpolation query points    - [nInterpPoints, nDim]
%
% Optional:
%  T   - delaunayn triangulation       - T = delaunayn(X);
%        Note, can be computed once and re-used for another data set V or for an
%        additional set of interpolation points XI in subsequent function calls.
%  method - 'linear','l' (default) OR 'nearest','n' - interpolation method
%  extrapval - extrapolation value, NaN (default)
%
% Output:
%  VI  - interpolated values           - [nInterpPoints, nData]
%  T   - delaunayn triangulation       - T = delaunayn(X);
%
% Comments:
%  Performs linear interpolation of a N-dimensional scattered data described
%  by the points X and values V. The interpolation is based on a N-dimensional
%  delaunayn triangulation. A barycentric interpolation scheme is employed for
%  all query points using tsearchn.
%  For slightly faster interpolation, choose nearest neighbour interpolation.
%  The nearst interpolation uses dsearchn instead of tsearchn.
%  Obs, 1-dimensional data is not supported, use interp1 instead.
%
%  Matlabs scatteredInterpolant class similarly allows for linear and nearest
%  neighbour scattered data interpolation. It is also significantly faster than
%  this function and have support for extrapolation. However, it can only handle
%  2D and 3D scatter data, whereas this function can handle any number of
%  dimensions.
%
%  Computationally, the main bottleneck in this function is the call to tsearchn
%  or dsearchn. Specifically because the underlying matlab mex-function tsrchnmx
%  needs transposed coordinate data.
%
% See also: delaunayn, tsearchn, dsearchn, scatteredInterpolant, interpn

%   Created by: Johan Winges
%   $Revision: 1.0$  $Date: 2016-01-11 14:00:00$


% Check for additional input:
narg = nargin;
ExtrapVal = NaN;
method_arg = 'linear';
if narg >= 4
   if isnumeric(varargin{end}) && isscalar(varargin{end}) && ischar(varargin{end-1})
      % User supplied an extrap val
      ExtrapVal = varargin{end};
      method_arg = varargin{end-1};
      narg = narg-2;
   elseif ischar(varargin{end})
      % User supplied a method
      method_arg = varargin{end};
      narg = narg-1;
   else
      method_arg = 'linear';
   end
end
if strncmpi(method_arg,'n',1)
   method_arg = 'n';
elseif strncmpi(method_arg,'l',1)
   method_arg = 'l';
end


if narg == 4
   % User supplied a triangulation T:
   T = varargin{1};
else
   T = [];
end

% Check dimension:
nDim = size(X, 2);
if nDim == 1 
   % Linear interpolation
   warning('scatteredInterpolant called with one-dimensional data, use interp1 instead')
   VI = interp1(X, V, XI, method_arg, ExtrapVal);
   return
end

% Make nd-triangulation if not supplied:
if isempty(T)
   T = delaunayn(X);
end

% Allocate output as extrapval:
sv  = size(V);
nXI = size(XI, 1);
VI  = ExtrapVal.*ones(cat(2,nXI, sv(2:end)));

if method_arg == 'l'
   % Find simplex and barycentric coordinates:
   [Tidx, BC] = tsearchn(X, T, XI);
   % Barycentric interpolation:
   for i = 1:nXI
      if ~isnan(Tidx(i))
         thisTidx = T(Tidx(i), :);
         thisV = V(thisTidx, :,:);
         thisBC = BC(i, :);
         thisVI = sum(bsxfun(@times, thisBC.', thisV), 1);
         VI(i, :) = thisVI(:);
      end
   end
elseif method_arg == 'n'
   % Find indices of closest point in X:
   k = dsearchn(X,T,XI,nan);
   VI(~isnan(k),:) = V(k(~isnan(k)),:);
end