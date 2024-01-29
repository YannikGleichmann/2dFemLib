classdef Mesh
  properties
    p                   % point matrix nDof x 2
    t                   % element/triangle matrix nTri x dim(V)
    e                   % edge matrix
    e_dirichlet         % edge matrix for diriclet bnd
    e_neumann           % edge matrix for neumann bnd
    e_src               % edge matrix for sommerfeld radiation cond
    e_abc               % edge matrix for absorbing boundary cond
    FunctionSpace       % function space i.e. P1, P1b, P2, P2b, ...
    nTriInit            % number of initial triangles (for coarsening); works only for P1 (or P1b) meshes!
    optional1
    optional2
    optional3
    optional4
    optional5
  end
  
  methods
  %% numbering
    function n = numberDof(obj)
      n = size(obj.p,1);
    end
    function n = numberTri(obj)
      n = size(obj.t,1);
    end
    function n = numberEdge(obj)
      n = size(obj.e,1);
    end
    function n = numberEdgeDirichlet(obj)
      n = size(obj.e_dirichlet,1);
    end
    function n = numberEdgeNeumann(obj)
      n = size(obj.e_neumann,1);
    end
    function n = numberEdgeSRC(obj)
      n = size(obj.e_src,1);
    end
    function n = numberEdgeABC(obj)
      n = size(obj.e_abc,1);
      end
  %% mesh object out of point and element/triangle matrix
      function obj = Mesh(p,t,nTriInit,FunctionSpace)
        if nargin > 0
          obj = Mesh;
          obj.p = p;
          obj.t = t;
          obj.e = obj.EdgeMatrix;
          if ~exist('nTriInit','var') || isempty('nTriInit')
            nTriInit = length(obj.t);
          end
          % determine FunctionSpace from size of t!
          if ~exist('FunctionSpace','var') || isempty('FunctionSpace')
            FunctionSpace = 'P1';
          end
          obj.nTriInit = nTriInit;
          obj.FunctionSpace = FunctionSpace;
        end
      end
  %% mesh quantities
      %% mesh element specs
      % computes the largest outer radius (R), inner radius (r), largest
      % side length (H), and smalles side length (h) of every triangle
      function [R,r,H,h,area] = Radii(obj)
        P1 = obj.p(obj.t(:,1), :);
        P3 = obj.p(obj.t(:,2), :);
        P2 = obj.p(obj.t(:,3), :);
        
        v1 = P2-P1; v2 = P3-P1;
        DetJ_K = v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2);
        AbsDetJ_K = abs(DetJ_K);
        
        a = P2-P3;
        b = P1-P3;
        c = P2-P1;
        len_a = sqrt(a(:,1).^2 + a(:,2).^2);
        len_b = sqrt(b(:,1).^2 + b(:,2).^2);
        len_c = sqrt(c(:,1).^2 + c(:,2).^2);
        
        R = len_a.*len_b.*len_c./AbsDetJ_K/2;
        if nargout <= 1; return; end
        r = AbsDetJ_K./(len_a + len_b + len_c);
        if nargout <= 2; return; end
        H = max([len_a,len_b,len_c],[],2);
        if nargout <= 3; return; end
        h = min([len_a,len_b,len_c],[],2);
        if nargout <= 4; return; end
        area = AbsDetJ_K(:);
      end
      %% boundary and interior nodes
      % markes the boundary (bnd) and interior (int) nodes
      function [bnd,int] = idxBoundaryInteriorEdges(obj)
        [r,b] = FunctionSpace2Dim(obj);
        if r == 1
          PlotMatrix = [1,2,3];
        elseif b
          r = r-1;
          if r == 2
            PlotMatrix = [1,4,7; 1,6,7; 3,5,6; 2,4,5; 4,5,7; 5,6,7];
          elseif r == 3
            PlotMatrix = [1,4,10; 1,9,10; 2,5,6; 3,7,8; 4,5,11; 4,10,11;
                          5,6,11; 6,7,11; 7,8,12; 7,11,12; 8,9,12;
                          9,10,12; 10,11,12];
          else; error('not implemented');
          end
        elseif ~b
          if r == 2
            PlotMatrix = [1,4,6; 2,4,5; 3,5,6; 4,5,6];
          elseif r == 3
            PlotMatrix = [1,9,10; 1,4,10; 2,5,6; 3,7,8; 4,5,6; 4,6,10;
                          6,7,10; 7,9,10; 7,8,9];
          else; error('not implemented');
          end
        else; error('FunctionSpace not recognized');
        end
        T = obj.t;
        P = obj.p;
        T = reshape(T(:,PlotMatrix')',3,[])';
        % own routine
%           P1 = P(T(:,1),:); P2 = P(T(:,2),:); P3 = P(T(:,3),:);
%           v1 = P2-P1; v2 = P3-P1;
%           idx = (v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2))<0;
%           T(idx,:) = T(idx,[2,1,3]);
%           E = [T(:,[1,2]); T(:,[2,3]); T(:,[3,1])];
%           re = E(:,[2,1]);
%           [~,pos] = intersect(E, re,'rows','stable');
%           bnd_e_idx = setdiff(1:size(E,1),pos);
%           E = E(bnd_e_idx,:);
%           bnd = unique([E(:,1) E(:,2)]);          % boundary nodes
%           int = setdiff(1:numberDof(obj),bnd);    % interior nodes
        % matlab function
        TR = triangulation(T(:,1:3),P);
        bnd = freeBoundary(TR); bnd = unique(bnd(:)); bnd = bnd(:);
        int = setdiff(1:numberDof(obj),bnd); int = int(:);
      end
      function e = EdgeMatrix(obj)
        e = [obj.t(:,[1,2]); obj.t(:,[2,3]); obj.t(:,[3,1])];
        re = e(:,[2,1]);
        [~,pos] = intersect(e, re, 'rows', 'stable');
        bnd_e_idx = setdiff(1:size(e,1),pos);
        e = e(bnd_e_idx,:);
      end
      %% function space specs
      % finds the dimension (r) of the function space; if b = 0 there are
      % no bubble functions, if b = 1 there are bubble functions
      function [r,b] = FunctionSpace2Dim(obj)
        if strcmp(obj.FunctionSpace,'P1')
          r = 1; if nargout<=1; return; end
          b = 0;
        elseif strcmp(obj.FunctionSpace,'P2')
          r = 2; if nargout<=1; return; end
          b = 0;
        elseif strcmp(obj.FunctionSpace,'P3')
          r = 3; if nargout<=1; return; end
          b = 0;
        elseif strcmp(obj.FunctionSpace,'P1b')
          r = 1; if nargout<=1; return; end
          b = 1;
        elseif strcmp(obj.FunctionSpace,'P2b')
          r = 3; if nargout<=1; return; end
          b = 1;
        elseif strcmp(obj.FunctionSpace,'P3b')
          r = 4; if nargout<=1; return; end
          b = 1;
        else
          error('not implemented');
        end
      end
      %% plot routine
      % a quick and dirty plot routine for Mesh; plots only the P1
      % equivalent mesh
      function fig = Plot(obj)
        fig = patch('Faces',obj.t(:,1:3),'Vertices',obj.p,...
               'FaceColor','none','EdgeColor','black','LineWidth',1);
      end
  end
end