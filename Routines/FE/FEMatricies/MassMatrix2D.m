function M = MassMatrix2D(mesh,m,lump,QuadratureOrder)
% this function computes the mass matrix
%       int_Omega m(x) phi_j(x) phi_i(x) dx
% 
%--------------------------------------------------------------------------
% Usage:
%           M = MassMatrx2D(mesh)
%           M = MassMatrx2D(mesh,m)
%           M = MassMatrx2D(mesh,m,lump)
%
% Input:
% * mesh        - mesh of Mesh class
% * m           - function of Function class or function handle
% * lump        - if lump = 0, mass matrix will NOT be lumped (only active
%                 when using bubble functions
%                 default: lump = 1;
% 
% Output:
% * M           - FE mass matrix (sparse) array nDof x nDof
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   21.10.2022

% folder = fileparts(which(mfilename));
% addpath(genpath(folder));

if exist('lump','var') == 0 || isempty(lump); lump = 1; end
% properties (dimension, bubblespace) of FunctionSpace and variation
[r,b] = FunctionSpace2Dim(mesh);
dim = 2*r;
if exist('m','var') == 0 || isempty(m)
    dim = dim+0;
elseif isa(m,'function_handle')
    dim = dim+2*r;
elseif isobject(m)
    rm = FunctionSpace2Dim(m.Mesh);
    dim = dim+rm;
end

% choose quadrature rule
if ~b
    tbl = Order2Quadrature(dim+1);
elseif ~lump
    tbl = Order2Quadrature(dim+1);
elseif b
    tbl = Order2Quadrature(r,b);
end

if exist('QuadratureOrder','var'); tbl = Order2Quadrature(QuadratureOrder); end

p = mesh.p;                 % point list
t = mesh.t;                 % element/triangle list
nDof = numberDof(mesh);     % number of points
nTri = numberTri(mesh);     % number of elements/triangle

P1 = p(t(:,1),:);           % corner points of each element/triangle
P2 = p(t(:,2),:);
P3 = p(t(:,3),:);

% quadratur
xi = tbl(1,:); eta = tbl(2,:); w = tbl(3,:);    % xi = x-node, eta = y-node, w = weight
nQuad = length(w);
loc_qnodes_x = ones(nTri,1)*xi;                                     % local
loc_qnodes_y = ones(nTri,1)*eta;                                    % quadrature nodes
glob_qnodes_x = P1(:,1)*(1-xi-eta) + P2(:,1)*xi + P3(:,1)*eta;      % global
glob_qnodes_y = P1(:,2)*(1-xi-eta) + P2(:,2)*xi + P3(:,2)*eta;      % quadrature nodes

% affine transformation F_K(xi) = J_K*xi + P1, where
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
if exist('m','var') == 0 || isempty(m)
    mQ = 1;
elseif isa(m,'function_handle')                 % a is in mesh.FunctionSpace (same FunctionSpace as the basis functions)
    mQ = m(glob_qnodes_x,glob_qnodes_y);        % so evaluate them at the quadrature nodes
elseif isa(m,'numeric')
    m = m(:);
    if size(mesh.p,1)~=size(m,1);error('wrong dimension of variation'); end
    m = Function(mesh,m);
    mQ = Interpolate(m,[glob_qnodes_x(:),glob_qnodes_y(:)]);
    mQ = reshape(mQ,nTri,nQuad);
elseif isobject(m)
    if isa(m.NodalValues,'function_handle')     % a is a function handle, so evaluate it at the grid points;
                                                % not at the Qnodes since it would be too exact; however, who would do that?
        m.NodalValues = m.NodalValues(m.Mesh.p(:,1), m.Mesh.p(:,2));
    end
    mQ = Interpolate(m,[glob_qnodes_x(:),glob_qnodes_y(:)]);
    mQ = reshape(mQ,nTri,nQuad);
end

% local assembling
M_loc = zeros(nTri,nBasis,nBasis);
for i = 1:nBasis
    for j = 1:nBasis
        M_loc(:, i, j) = AbsDetJ_K.*((mQ.*(N{i}.*N{j}))*w(:));
    end
end

% global mass matrix
[elementindex,col,row] = ndgrid(1:nTri, (1:nBasis)-1, (1:nBasis)-1);   % Zeilen-/Spaltenindizes
i = t(elementindex + row*nTri);                 % Zeilenindex
j = t(elementindex + col*nTri);                 % Spaltenindex
M = sparse(i(:), j(:), M_loc(:), nDof, nDof);   % Globale Steifigkeitsmatrix
if b && lump
    M = sum(M,2);                   % in M there could be too small values
    M = spdiags(M,0,nDof,nDof);     % so "get rid" of them
end
M = (M+M')/2;                       % force M to be symmetric (due to computation errors)