function mesh = MeshFromEdges(x,y,orientation,varargin)
% This function generates a mesh given from edge points x_i, y_i for P1
% elements. For higher order meshes use the genMesh function after
% this function (e.g. with nDofMax set to 0).
%
% Usage:
%       mesh = MeshFromEdges(x,y,orentation)
%       mesh = MeshFromEdges(x,y,orentation,varargin)
%
% -------------------------------------------------------------------------
% Input:
%  * x,y             - edege points x_i, y_i that define the different
%                      line segments from points x_i,y_i to x_i+1,y_i+1
%  * orientation     - 'clockwise' or 'counterclockwise':
%                       going from the i-th to the i+1-th point clockwise
%                       or counterclockwise
%  * varargin        - the same varargins as the Matlab function
%                      initmesh takes
% 
% -------------------------------------------------------------------------
% Output:
%  * Mesh object
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   18.04.2021
if isempty(varargin) || ~exist('varargin','var')
  varargin = {'Hmax',0.1};
end

x = x(:);
y = y(:);
nLineSegments = length(x);
% line segment
seg = 2*ones(nLineSegments,1);
% orientation
switch orientation
  case 'clockwise'
    or = ones(nLineSegments,1)*[0,1];
  case 'counterclockwise'
    or = ones(nLineSegments,1)*[1,0];
  otherwise
    warning('wrong orientation');
end
% geometry matrix
xPart = [x, [x(2:end); x(1)]];
yPart = [y, [y(2:end); y(1)]];
g = [seg,xPart,yPart,or]';
% construct mesh
[p,~,t] = initmesh(g,varargin{:});
% convert to mesh class
p = p';
t = t(1:3,:)';
mesh = Mesh(p,t);
