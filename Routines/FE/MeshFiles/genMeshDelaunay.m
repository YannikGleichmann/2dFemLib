function mesh = genMeshDelaunay(xMin,xMax,Nx,yMin,yMax,Ny,FunctionSpace,...
                                varargin)
% This function generates a 'P1'='P1b', 'P2', 'P2b', 'P3', or 'P3b' finite
% element mesh using Matlab's delaunay function for a rectangular domain
% [xMin,xMax] x [yMin,yMax] with Nx and Ny points in x and y direction,
% respectively, where the verticies lie on an Cartesian grid.
% 
%--------------------------------------------------------------------------
% References:
% [1] G. Cohen et. al. | Higher Order Triangular Finite Elements with
%                        Mass Lumping for the Wave Equation
% [2] W. A. Mulder | Higher-Order Mass-Lumped Finite Elements
%                    for the Wave Equation
% 
%--------------------------------------------------------------------------
% Usage:
%       mesh = genMeshDelaunay()
%       mesh = genMeshDelaunay(xMin,xMax,Nx,yMin,yMax,Ny,FunctionSpace)
%       and more, see tutorial
% 
% -------------------------------------------------------------------------
% Input:
% * xMin, xMax, yMin, yMax  - boundary of rectengular domain
% * Nx, Ny                  - number of grid points in each direction
% * FunctionSpace           - 'P1', 'P1b', 'P2', 'P2b', 'P3', 'P3b'
%                             suffix b: bubble function for mass lumping
%                             see [1,2]
% * varargin                - points given as x and y or a p matrix that
%                             should be enforced to be in the mesh
% 
%--------------------------------------------------------------------------
% Output:
% * Mesh object
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   25.10.2023

if exist('xMin','var') == 0 || isempty(xMin); xMin = 0; end
if exist('xMax','var') == 0 || isempty(xMax); xMax = 1; end
if exist('Nx','var') == 0 || isempty(Nx); Nx = 100; end
if exist('yMin','var') == 0 || isempty(yMin); yMin = 0; end
if exist('yMax','var') == 0 || isempty(yMax); yMax = 1; end
if exist('Ny','var') == 0 || isempty(Ny); Ny = 100; end
if exist('FunctionSpace','var') == 0 || isempty(FunctionSpace)
    FunctionSpace = 'P1';
end
if ~ischar(FunctionSpace); error('FunctionSpace must be a string'); end
numVarArgs = length(varargin);
if numVarArgs == 1
    Value = varargin{1};
    if size(Value,1) == 2
        Value = Value';
        x = Value(:,1); y = Value(:,2);
    elseif size(Value,2) == 2
        x = Value(:,1); y = Value(:,2);
    elseif isempty(Value)
        x = [];
        y = [];
    else
        error('point matrix should be an nx2 or 2xn array');
    end
elseif numVarArgs == 2
    x = varargin{1}; y = varargin{2};
    x = x(:); y = y(:);
elseif numVarArgs == 0
  x = [];
  y = [];
else
  error('Wrong input of arguments')
end

% regular mesh whos veticies lie on an Cartesian grid
[X,Y] = ndgrid(linspace(xMin,xMax,Nx),linspace(yMin,yMax,Ny));
p = [X(:),Y(:);x,y];
p = unique(p,'rows','legacy');
t = delaunay(p);

% construct higher order meshes
switch FunctionSpace
    case 'P1'
    case 'P1b'
    case 'P2'
        barycentric_xyz = PrNewtonCotesNodes(2);
        [p,t] = P1toPr(p,t,2,barycentric_xyz);
    case 'P2b'
        barycentric_xyz = PrBubbleNodes(2);
        [p,t] = P1toPr(p,t,2,barycentric_xyz);
    case 'P3'
        barycentric_xyz = PrNewtonCotesNodes(3);
        [p,t] = P1toPr(p,t,3,barycentric_xyz);
    case 'P3b'
        barycentric_xyz = PrBubbleNodes(3);
        [p,t] = P1toPr(p,t,3,barycentric_xyz);
    otherwise
        error('Function space %s not implemented.',FunctionSpace);
end

% construct edge matrix
e = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];
re = e(:,[2,1]);
% check, if the edge has a reversed counterpart; if so it is an interior edge
[~,pos] = intersect(e,re,'rows','stable');
bnd_e_idx = setdiff(1:size(e,1),pos);
e = e(bnd_e_idx,:);

% mesh object
mesh = Mesh(p,t,size(t,1),FunctionSpace);
mesh.e = e;