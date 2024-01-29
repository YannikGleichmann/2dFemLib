function L = LoadVector2D(mesh,l)
% this function computes the load vector
%       int_Omega l(x)*phi_i(x) dx
% 
%--------------------------------------------------------------------------
% Usage:
%           L = LoadVector2D(mesh)
%           L = LoadVector2D(mesh,l)
%
% Input:
% * mesh        - mesh of Mesh class
% * l           - function of Function class or function handle
% 
% 
% Output:
% * L           - FE load vector (sparse) array nDof x 1
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   21.10.2022

% properties (dimension, bubblespace) of FunctionSpace and variation
[r,b] = FunctionSpace2Dim(mesh);
dim = r;
if exist('l','var') == 0 || isempty(l)
  dim = dim+0;
elseif isa(l,'function_handle')
  dim = dim+2*r;
elseif isobject(l)
  la = FunctionSpace2Dim(l.Mesh);
  dim = dim+la;
end

% choose quadrature rule; here we do not need a specific quadrature rule
% for the bubble functions
tbl = Order2Quadrature(dim+1);

p = mesh.p;                 % point list
t = mesh.t;                 % element/triangle list
nDof = numberDof(mesh);     % number of points
nTri = numberTri(mesh);     % number of elements/triangles

P1 = p(t(:,1),:);           % corner points of each element/triangle
P2 = p(t(:,2),:);
P3 = p(t(:,3),:);

% quadratur
xi = tbl(1,:); eta = tbl(2,:); w = tbl(3,:); 
nQuad = length(w);
loc_qnodes_x = ones(nTri,1)*xi;
loc_qnodes_y = ones(nTri,1)*eta;
glob_qnodes_x = P1(:,1)*(1-xi-eta) + P2(:,1)*xi + P3(:,1)*eta;
glob_qnodes_y = P1(:,2)*(1-xi-eta) + P2(:,2)*xi + P3(:,2)*eta;

% affine transformation F_K(z) = J_K*z + P1, where
%                       J_K = [P2 - P1, P3 - P1] = [v1, v2]
%                       F_K' = J_K
v1 = P2-P1; v2 = P3-P1;
DetJ_K = v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2);
AbsDetJ_K = abs(DetJ_K);

% basis functions; recall that r = dim(P). So r = dim(P1b) = 1 but
% r = dim(Pib) = i+1; BubbleBasis2D take i as input argument.
if b && (r==1)
  Basis = BubbleBasis2D(r);
elseif b
  Basis = BubbleBasis2D(r-1);
elseif ~b
  Basis = LagrangeBasis2D(r);
end
N = Basis(loc_qnodes_x,loc_qnodes_y);
nBasis = length(N);

% variation
if exist('l','var') == 0 || isempty(l)
  lQ = 1;
elseif isa(l,'function_handle')
  lQ = l(glob_qnodes_x,glob_qnodes_y);
elseif isa(l,'numeric')
  l = l(:);
  if size(mesh.p,1)~=size(l,1);error('wrong dimension of variation'); end
  l = Function(mesh,l);
  lQ = Interpolate(l,[glob_qnodes_x(:),glob_qnodes_y(:)]);
  lQ = reshape(lQ,nTri,nQuad);
elseif isobject(l)
  if isa(l.NodalValues,'function_handle')
    l.NodalValues = l.NodalValues(l.Mesh.p(:,1), l.Mesh.p(:,2));
  end
  lQ = Interpolate(l,[glob_qnodes_x(:),glob_qnodes_y(:)]);
  lQ = reshape(lQ,nTri,nQuad);
end

l_loc = zeros(nTri,nBasis);
for i = 1:nBasis
    l_loc(:,i) = AbsDetJ_K.*(lQ.*N{i})*w(:);
end
L = full(sparse(t,1,l_loc,nDof,1));