function basis = PrLagrangeBasis2D(r,barycentric)
%% nodal nodes
if exist('barycentric','var') == 0 || isempty(barycentric)
  barycentric = PrNewtonCotesNodes(r);
end
nodal = barycentric(:,2:3);

%% power
P = repmat((0:r)',1,r+1);
pow_y = triu(P +1);
pow_x = triu(P'+2)-pow_y;
pow = ([nonzeros(pow_x),nonzeros(pow_y)]-1);
[X,DegX] = ndgrid(nodal(:,1),pow(:,1));
[Y,DegY] = ndgrid(nodal(:,2),pow(:,2));

%% monom
A = (X.^DegX) .* (Y.^DegY);

%% coefficient
coefficient = inv(A);
coefficient = round(coefficient,...
                    ceil(14-log(max(abs(coefficient(:))))/log(10)));

px  = pow(:,1);
py  = pow(:,2);

% power order encode
pow_map = inversemapping(pow*[1;r+1]+1);
n = size(pow,1);

% Dx
d = [px-1,py]; d(px==0,:) = 0;
Dx = sparse(n,n); Dx(sub2ind([n,n],pow_map(d*[1;r+1]+1),(1:n)')) = px;

% Dy
d = [px,py-1]; d(py==0,:) = 0;
Dy = sparse(n,n); Dy(sub2ind([n,n],pow_map(d*[1;r+1]+1),(1:n)')) = py;

monomstr = sprintf('(x.^%d).*(y.^%d),',pow');
monom = [];
eval(['monom = @(x,y) [' monomstr '];']);
basis = @Lagrange2D;

  function [f,fx,fy,fxx,fyy,fxy] = Lagrange2D(x)
    m = monom(x(:,1),x(:,2));
    f = m*coefficient;
    
    if nargout >=2
      coefficient_x = Dx*coefficient;
      coefficient_y = Dy*coefficient;
      fx = m*coefficient_x;
      fy = m*coefficient_y;
    end
    if nargout >=4
      coefficient_xx = Dx*coefficient_x;
      coefficient_xy = Dx*coefficient_y;
      coefficient_yy = Dy*coefficient_y;
      
      fxx = m*coefficient_xx;
      fxy = m*coefficient_xy;
      fyy = m*coefficient_yy;
    end
  end

end