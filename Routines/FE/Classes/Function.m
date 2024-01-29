classdef Function < handle
  properties
    NodalValues     % nodal values on verticies
    Mesh            % mesh of Mesh class
    optional1
    optional2
    optional3
    optional4
    optional5
  end
  
  methods
    %% plot function
    % quick-and-dirty function plot; plots only the P1 equivalent
    % functions, that is only the function on the P1/edge nodes
    function fig = Plot(obj)
      varargin = {'facecolor','interp','edgecolor','none'};
      if any(size(obj.NodalValues) == 1)
        fig = trisurf(obj.Mesh.t(:,1:3),...
                obj.Mesh.p(:,1),obj.Mesh.p(:,2),...
                obj.NodalValues,varargin{:});
      else
        [~,MaxIt] = size(obj.NodalValues);
        F = obj.NodalValues;
        mi = min(F(:)); ma = max(F(:));
        for k = 1:MaxIt
          fig = trisurf(obj.Mesh.t(:,1:3),...
                  obj.Mesh.p(:,1),obj.Mesh.p(:,2),...
                  F(:,k),varargin{:});
          zlim([mi,ma]); drawnow;
        end
      end
    end
    %% make function object out of Mesh and NodalValues
    function obj = Function(mesh,NodalValues)
      if nargin > 0
        if ~isa(mesh,'Mesh')
          error('enter mesh object');
        end
        obj.Mesh = mesh;
        if isa(NodalValues,'numeric')
          if size(NodalValues,1) ~= size(mesh.p,1)
            warning('NodalValues not consistent with mesh');
            obj.NodalValues = NodalValues;
          else
            obj.NodalValues = NodalValues;
          end
        elseif isa(NodalValues,'function_handle')
          warning('NodalValues is a function handle')
          obj.NodalValues = NodalValues;
        else
          error('wrong data type for NodalValues');
        end
      end
    end
    %% arithmetic operations
    function obj = plus(u,v)
      if ~isa(u,'Function') || ~isa(v,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      if u.Mesh.p ~= v.Mesh.p
        error('mesh mismatch');
      end
      obj = Function(u.Mesh,u.NodalValues+v.NodalValues);
    end
    function obj = minus(u,v)
      if ~isa(u,'Function') || ~isa(v,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      if u.Mesh.p ~= v.Mesh.p
        error('mesh mismatch');
      end
      obj = Function(u.Mesh,u.NodalValues-v.NodalValues);
    end
    function obj = times(u,v)
      if ~isa(u,'Function') || ~isa(v,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      if u.Mesh.p ~= v.Mesh.p
        error('mesh mismatch');
      end
      obj = Function(u.Mesh,u.NodalValues.*v.NodalValues);
    end
    function obj = divide(u,v)
      if ~isa(u,'Function') || ~isa(v,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      if u.Mesh.p ~= v.Mesh.p
        error('mesh mismatch');
      end
      if any(v.NodalValues(:) == 0)
        error('no division by 0 possible');
      end
      obj = Function(u.Mesh,u.NodalValues./v.NodalValues);
    end
    function obj = power(u,n)
      if ~isa(u,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      if ~isa(n,'numeric')
        error('numeric value for exponent expected');
      end
      obj = Function(u.Mesh,u.NodalValues.^n);
    end
    function obj = sqrt(u)
      if ~isa(u,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      obj = Function(u.Mesh,sqrt(u.NodalValues));
    end
    %% error and norms
    function n = L2norm(u)
      if ~isa(u,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      M = MassMatrix2D(u.Mesh);
      n = sqrt(sum(u.NodalValues.*(M*u.NodalValues),1));
    end
    function n = semiH1norm(u)
      if ~isa(u,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      A = StiffnessMatrix2D(u.Mesh);
      n = sqrt(sum(u.NodalValues.*(A*u.NodalValues),1));
    end
    function n = L2error(u,v)
      if ~isa(u,'Function') || ~isa(v,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      if u.Mesh.p ~= v.Mesh.p
        error('mesh mismatch');
      end
      M = MassMatrix2D(u.Mesh);
      n = sqrt((u.NodalValues-v.NodalValues)'*...
                M*(u.NodalValues-v.NodalValues));
    end
    function [H1_loc,H1_locScaled] = semiH1normElement(u)
      if ~isa(u,'Function')
        error('wrong input data type; Only Function objects are allowed');
      end
      if ~any(size(u.NodalValues) == 1)
        error('wrong dimension of nodal values');
      end
      mesh = u.Mesh; U = u.NodalValues;
      [r,b] = FunctionSpace2Dim(mesh); dim = 2*(r-1);
      tbl = Order2Quadrature(dim+1);
      p = mesh.p; t = mesh.t;
      nTri = numberTri(mesh);
      P1 = p(t(:,1),:); P2 = p(t(:,2),:); P3 = p(t(:,3),:);
      xi = tbl(1,:); eta = tbl(2,:); w = tbl(3,:);
      loc_qnodes_x = ones(nTri,1)*xi;
      loc_qnodes_y = ones(nTri,1)*eta;
      v1 = P2-P1; v2 = P3-P1;
      DetJ_K = v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2);
      AbsDetJ_K = abs(DetJ_K);
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
        dNdx{i} = dNdxi{i}.* (v2(:,2)./DetJ_K) + dNdeta{i}.*...
                  (-v1(:,2)./DetJ_K);
        dNdy{i} = dNdxi{i}.* (-v2(:,1)./DetJ_K) + dNdeta{i}.*...
                  (v1(:,1)./DetJ_K);
      end
      A_loc = zeros(nTri, nBasis, nBasis);
      for i = 1:nBasis
        for j = 1:nBasis
          A_loc(:, i, j) = AbsDetJ_K.*(((dNdx{i}.*dNdx{j} ...
                                         + dNdy{i}.*dNdy{j}) )*w(:));
        end
      end
      U_loc = zeros(nTri,nBasis);
      for i = 1:nBasis
        U_loc(:,i) = U(t(:,i));
      end
      U_loc = U_loc';
      H1_loc = zeros(nTri,1);
      for k = 1:nTri
        H1_loc(k) = U_loc(:,k)'*squeeze(A_loc(k,:,:))*U_loc(:,k);
      end
      H1_loc = sqrt(H1_loc);
      if nargout > 1; H1_locScaled = H1_loc.*AbsDetJ_K; end
    end
  end
end