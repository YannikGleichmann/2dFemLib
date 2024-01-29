function A = StiffnessMatrix2D(mesh,a)
% this function computes the stiffness matrix
%       int_Omega a(x) Dphi_j(x) x Dphi_i(x) dx
% 
%--------------------------------------------------------------------------
% Usage:
%           A = StiffnessMatrix2D(mesh)
%           A = StiffnessMatrix2D(mesh,a)
%
% Input:
% * mesh        - mesh of Mesh class
% * a           - function of Function class or function handle
% 
% 
% Output:
% * A           - FE stiffness matrix (sparse) array nDof x nDof
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   21.10.2022

% properties (dimension, bubblespace) of FunctionSpace and variation
[r,b] = FunctionSpace2Dim(mesh);
dim = 2*(r-1);
if exist('a','var') == 0 || isempty(a)
    dim = dim+0;
elseif isa(a,'function_handle')
    dim = dim+2*r;
elseif isobject(a)
    ra = FunctionSpace2Dim(a.Mesh);
    dim = dim+ra;
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
[N,dNdxi,dNdeta] = Basis(loc_qnodes_x,loc_qnodes_y);
nBasis = length(N);
dNdx = cell(size(dNdxi));
dNdy = cell(size(dNdeta));

for i = 1:nBasis
  dNdx{i} = dNdxi{i}.* (v2(:,2)./DetJ_K) + dNdeta{i}.* (-v1(:,2)./DetJ_K);
  dNdy{i} = dNdxi{i}.* (-v2(:,1)./DetJ_K) + dNdeta{i}.* (v1(:,1)./DetJ_K);
end

% variation
if exist('a','var') == 0 || isempty(a)
  aQ = 1;
elseif isa(a,'function_handle')
  aQ = a(glob_qnodes_x,glob_qnodes_y);
elseif isa(a,'numeric')
  a = a(:);
  if size(mesh.p,1)~=size(a,1);error('wrong dimension of variation'); end
  a = Function(mesh,a);
  aQ = Interpolate(a,[glob_qnodes_x(:),glob_qnodes_y(:)]);
  aQ = reshape(aQ,nTri,nQuad);
elseif isobject(a)
  if isa(a.NodalValues,'function_handle')
    a.NodalValues = a.NodalValues(a.Mesh.p(:,1), a.Mesh.p(:,2));
  end
  aQ = Interpolate(a,[glob_qnodes_x(:),glob_qnodes_y(:)]);
  aQ = reshape(aQ,nTri,nQuad);
end

% local assembling
A_loc = zeros(nTri, nBasis, nBasis);
for i = 1:nBasis
  for j = 1:nBasis
    A_loc(:, i, j) = AbsDetJ_K.*(( aQ.* ...
                     (dNdx{i}.*dNdx{j} + dNdy{i}.*dNdy{j}) )*w(:));
  end
end

% global assembling
[elementindex, col, row] = ndgrid(1:nTri, (1:nBasis) - 1, (1:nBasis) - 1);
i = t(elementindex + row*nTri);
j = t(elementindex + col*nTri);
A = sparse(i(:), j(:), A_loc(:), nDof, nDof);
A = (A+A')/2;