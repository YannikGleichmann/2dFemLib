function mesh = genMesh(varargin)
% This function generates a 'P1'='P1b', 'P2', 'P2b', 'P3', or 'P3b' finite
% element mesh
% 
%--------------------------------------------------------------------------
% References:
% [1] G. Cohen et. al. | Higher Order Triangular Finite Elements with
%                        Mass Lumping for the Wave Equation
% [2] W. A. Mulder | Higher-Order Mass-Lumped Finite Elements
%                    for the Wave Equation
% [3] Stefan A. Funken, Anja Schmidt | Adaptive Mesh Refinement in 2D -
%                                      An Efficient Implementation in
%                                      Matlab
%
% [4] https://github.com/dengwirda/mesh2d
% 
% -------------------------------------------------------------------------
% Usage:
%       mesh = genMesh()
%       mesh = genMesh(Parameter,Value)
%       mesh = genMesh(Parameter1,Value1,Parameter2,Value2,...)
% 
%           Parameter, Value: 'PointMatrix', p
%                             'ElementMatrix', t
%                             'FunctionSpace', 'P' (string)
%                             'nRef', integer >= 0
%                             'hmax': float > 0
% 
% -------------------------------------------------------------------------
% Input:
% * p               - initial point matrix   (nDof x 2) \ has to be for P1
% * t               - initial element matrix (nTri x 3) / elements
% * FunctionSpace   - 'P1', 'P1b', 'P2', 'P2b', 'P3', 'P3b' elements;
%                     suffix b: bubble function for mass lumping, see [1]
% * hmax            - maximal smallest side length
% * nRef            - number of refinements from original mesh
% 
%--------------------------------------------------------------------------
% Output:
% * Mesh object
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   25.10.2022

numVarArgs = length(varargin);
if mod(numVarArgs,2) && (numVarArgs > 0)
  error('Check input arguments.');
end
cnt = 0;
% read input arguments
for ii = 1:2:numVarArgs
  Param = varargin{ii}; Value = varargin{ii+1};
  if strcmp(Param,'PointMatrix')
    p = Value;
  elseif strcmp(Param,'ElementMatrix')
    t = Value;
  elseif strcmp(Param,'FunctionSpace')
    FunctionSpace = Value;
  elseif strcmp(Param,'nRef')
    nRef = Value; cnt = cnt + 1;
  elseif strcmp(Param,'hmax')
    hmax = Value; cnt = cnt + 1;
  end
end

if cnt == 0
  hmax = 0.1;
  fprintf('No nRef or hmax specified. \n Use hmax = 0.1.\n');
elseif cnt >= 2
  error('specify only one parameter: nRef or hmax');
end

% check if point and element matrix exist; if not construct mesh on [0,1]^2
if ~exist('p','var') && ~exist('t','var')
  fprintf('constructing mesh on [0,1] x [0,1].\n');
  p = [0,0; 1,0; 1,1; 0,1];   % create initial mesh on
  t = [1,2,4; 2,3,4];         % [0,1] x [0,1]
end

% check FunctionSpace; if not specifield use P1
if ~exist('FunctionSpace','var')
  FunctionSpace = 'P1';
  fprintf('No function space specified. P1 will be used.\n');
else
  if ~ischar(FunctionSpace); error('FunctionSpace must be a string'); end
end
 
% construct P1 mesh
if exist('hmax','var')
  hdata.hmax = hmax;
  warning off;
  [p,t] = mesh2d(p,[],hdata);
  warning on;
elseif exist('nRef','var')
  if ~exist('p','var') || ~exist('t','var')
    error('Specify both, point and element matrix.');
  end
  for n = 1:nRef
%       [p,t] = TrefineRGB(p,t,1:size(t,1));
    [p,t] = refine(p,t,true(size(t,1),1));
  end
end
    
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
% check, if the edge has a reversed counterpart;
% if so it is an interior edge
[~,pos] = intersect(e,re,'rows','stable');
bnd_e_idx = setdiff(1:size(e,1),pos);
e = e(bnd_e_idx,:);

% mesh object
mesh = Mesh(p,t,size(t,1),FunctionSpace);
mesh.e = e;