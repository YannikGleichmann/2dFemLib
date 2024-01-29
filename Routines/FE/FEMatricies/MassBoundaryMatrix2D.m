function S = MassBoundaryMatrix2D(mesh,s,lump,E)
% this function computes the mass matrix on the boundary Gamma
%       int_Gamma m(x) phi_j(x) phi_i(x) dx
% 
%--------------------------------------------------------------------------
% Usage:
%           S = MassMatrx2D(mesh)
%           S = MassMatrx2D(mesh,s)
%           S = MassMatrx2D(mesh,s,lump)
%
% Input:
% * mesh        - mesh of Mesh class
% * s           - function of Function class or function handle
% * lump        - if lump = 0, mass matrix will NOT be lumped (only active
%                 when using bubble functions
%                 default: lump = 1;
% 
% Output:
% * S           - FE matrix (sparse) array nDof x nDof
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   21.10.2022

if exist('lump','var') == 0 || isempty(lump); lump = 1; end
% properties (dimension, bubblespace) of FunctionSpace and variation
[r,b] = FunctionSpace2Dim(mesh);
dim = 2*r;
if exist('s','var') == 0 || isempty(s)
    dim = dim+0;
elseif isa(s,'function_handle')
    dim = dim+r;
elseif isobject(s)
    rs = FunctionSpace2Dim(s.Mesh);
    dim = dim+rs;
end

% choose quadrature rule
if ~b
    tbl = Order2Quadrature1D(dim+1);
elseif ~lump
    tbl = Order2Quadrature1D(dim+1);
elseif b
    tbl = Order2Quadrature1D(r,b);
end

p = mesh.p;                 % point list
t = mesh.t;                 % element/triangle list
e = mesh.e;                 % triangle list
if exist('E','var') && ~isempty(E); e = E; end

nDof = numberDof(mesh);     % number of points
nTri = numberTri(mesh);     % number of elements/triangles
nEdges = length(e);  % number of edges

t12 = t(:,1:2);
t23 = t(:,2:3);
t31 = t(:,[3,1]);

[idx_12, tri12] = ismember(e, t12, 'rows', 'legacy');
[idx_23, tri23] = ismember(e, t23, 'rows', 'legacy');
[idx_31, tri31] = ismember(e, t31, 'rows', 'legacy');
idx_elements = [tri12; tri23; tri31];               % the i-th entry maps the i-th
idx_elements(~[idx_12; idx_23; idx_31]) = [];       % edge to its corresponding triangle
idx_sides = [1*ones(nnz(idx_12), 1); ...            % ---------------------
             2*ones(nnz(idx_23), 1); ...            % in {1, 2, 3}; idx_sides(i) marks the side (first, second, third)
             3*ones(nnz(idx_31), 1)];               % in which the i-th edge is on its corresponding triangle
lenEdges = sqrt(sum(( p(e(:, 1), :) - p(e(:, 2), :) ).^2, 2));

% quadratur, parametrization
xi = tbl(1,:); xi = xi(:);
w = tbl(2,:); w = w(:);

nQuad = length(w);
% 3   unit triangle; parametrization g(t): (0, 1) -> (x, y)
% |\        from (1) to (2): g(t) = (t, 0) = (X1(1), Y1(1))
% | \       from (2) to (3):      = (1-t, t)
% |  \      from (3) to (1):      = (0, 1-t) = (X1(3), Y1(3))
% 1---2
O = 0*xi;
X1 = [xi, 1-xi, O];
Y1 = [O, xi, 1-xi];
x = X1(:, idx_sides)';
y = Y1(:, idx_sides)';

% basis functions; recall that r = dim(P). So r = dim(P1b) = 1 but
% r = dim(Pib) = i+1; BubbleBasis2D take i as input argument.
if b && (r==1)
    Basis = BubbleBasis2D(r);
elseif b
    Basis = BubbleBasis2D(r-1);
elseif ~b
    Basis = LagrangeBasis2D(r);
end
phi = Basis(x,y);
nBasis = length(phi);


% variation
if exist('s','var') == 0 || isempty(s)
    sQ = 1;
elseif isa(s,'function_handle')                 % a is in mesh.FunctionSpace (same FunctionSpace as the basis functions)
    sQ = s(x,y);        % so evaluate them at the quadrature nodes
elseif isa(s,'numeric')
    s = s(:);
    if size(mesh.p,1)~=size(s,1);error('wrong dimension of variation');end
    s = Function(mesh,s);
    sQ = Interpolate(s,[x(:),y(:)]);
    sQ = reshape(sQ,nEdges,nQuad);
elseif isobject(s)
    if isa(s.NodalValues,'function_handle')     % a is a function handle, so evaluate it at the grid points;
                                                % not at the Qnodes since it would be too exact; however, who would do that?
        s.NodalValues = s.NodalValues(s.Mesh.p(:,1), s.Mesh.p(:,2));
    end
    sQ = Interpolate(s,[x(:),y(:)]);
    sQ = reshape(sQ,nEdges,nQuad);
end

% local assembling
S_loc = zeros(nEdges,nBasis,nBasis);
for i = 1:nBasis
  for j = 1:nBasis
    S_loc(:,i,j) = lenEdges.*((sQ.*(phi{j}.*phi{i}))*w(:));
  end
end

% global assembling
[EDindex, I, J] =  ndgrid(idx_elements, (1:nBasis)-1, (1:nBasis)-1);
EL_i_idx = t(EDindex + I*nTri);
EL_j_idx = t(EDindex + J*nTri);
S = sparse(EL_i_idx(:), EL_j_idx(:), S_loc(:), nDof, nDof);

if b && lump
    S = sum(S,2);                   % in S there could be values < 1e-30
    S = spdiags(S,0,nDof,nDof);     % so "get rid" of them
end
S = (S+S')/2;                       % force S to be symmetric (due to computation errors)