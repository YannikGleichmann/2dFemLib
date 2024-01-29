function fInter = Interpolate(f,varargin)
% this functions interpolates the FE function f to query points
% 
%--------------------------------------------------------------------------
% References regarding triangles:
% https://ch.mathworks.com/help/matlab/ref/triangulation.html
% 
%--------------------------------------------------------------------------
% Usage:
%       fInter = Interpolate(f,qp)
%       fInter = Interpolate(f,qp_x,qp_y);
% 
% Input:
% * f               - (interpoland) object of Function class
% * qp              - n x 2 array with the query points
%   or qp_x,qp_y    - x and y value of the query points
% 
% Output:
% * fInter          - f interpolated at the query points
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   13.06.2020

% folder = fileparts(which(mfilename));
% addpath(genpath(folder));
if ~isobject(f); error('function must be a Function class'); end
numVarArgs = length(varargin);
if numVarArgs == 0
    error('specify query points');
elseif numVarArgs == 1
    Value = varargin{1};
    if size(Value,1) == 2
        Value = Value';
        x = Value(:,1); y = Value(:,2);
    elseif size(Value,2) == 2
        x = Value(:,1); y = Value(:,2);
    else
        error('point matrix should be an nx2 or 2xn array');
    end
elseif numVarArgs == 2
    x = varargin{1}; y = varargin{2};
    x = x(:); y = y(:);
end
% deterime FunctionSpace for f and choose correct asis
[r,b] = FunctionSpace2Dim(f.Mesh);
if b && (r==1)
    Basis = BubbleBasis2D(r);
elseif b
    Basis = BubbleBasis2D(r-1);
elseif ~b
    Basis = LagrangeBasis2D(r);
end
% cosmetic for readability
fNodes = f.NodalValues(:);
t = f.Mesh.t;
p = f.Mesh.p;
% check in which triangle the points in p lie
% ID(i) = k says that the i-th point in p lies in the k-th triangle/element
warning off; TR = triangulation(t(:,1:3),p);
ID = pointLocation(TR,x,y); warning on;
% compute the local counterpart of the points
% i.e. F_K(x_loc) = [x;y], where F_K(xi) = J_K*xi + P1
P1 = p(t(ID,1),:);      % corner points of the triangle to
P2 = p(t(ID,2),:);      % which the i-th point belongs
P3 = p(t(ID,3),:);
% compute F_K, then compute xLoc = J_K^-1(x-P1)
v1 = P2-P1; v2 = P3-P1;
DetJ_K = v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2);
P = [x,y]-P1;
xLoc = [v2(:,2).*P(:,1) - v2(:,1).*P(:,2), ...
         -v1(:,2).*P(:,1)+v1(:,1).*P(:,2)]./DetJ_K;
% evaluate Basis functions at the local points
N = Basis(xLoc(:,1),xLoc(:,2));
nBasis = length(N);
fInter = zeros(size(x));
for i = 1:nBasis
    fInter = fInter + fNodes(t(ID,i)).*N{i};
end