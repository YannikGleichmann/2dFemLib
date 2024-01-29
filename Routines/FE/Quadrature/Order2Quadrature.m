function tbl = Order2Quadrature(r,b)
% this function return a quadrature rule of order r as a table tbl where
% row 1-2 corresponds to the x and y value, respectively and row 3 are the
% weights to the corresponding quadrature points.
% 
% -------------------------------------------------------------------------
% Usage:
%       tbl = Order2Quadrature(r)
%       tbl = Order2Quadrature(r,[]) ( = Order2Quadrature(r) )
%       tbl = Order2Quadrature(r,b)
% 
% -------------------------------------------------------------------------
% Input:
%   * r     - order of Gaussian quadrature rule
%   * b     - =1 if quadrature rule for mass lumping is desired
% 
% -------------------------------------------------------------------------
% Output:
%   * tbl = xi(1),   xi(2), ..., xi(nQuad)      <- x values
%           eta(1), eta(2), ..., eta(nQuad)     <- y values
%           w(1),     w(2), ..., w(nQuad)       <- weights
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   14.06.2020

if exist('b','var') == 0 || isempty(b); b = 0; end
triangle = QuadratureCoefficientTriangle();

if b
  switch r
    case 1
      tbl = triangle.Bubble3;
    case 2
      error('no bubble space with dimension %d',r');
    case 3
      tbl = triangle.Bubble7;
    case 4
      tbl = triangle.Bubble12;
    otherwise
      error('not implemented');
  end
else
  switch r
    case 1
      tbl = triangle.Gauss1;
    case 2
      tbl = triangle.Gauss3;
    case 3
      tbl = triangle.Gauss6;
    case 4
      tbl = triangle.Gauss6;
    case 5
      tbl = triangle.Gauss7;
    case 6
      tbl = triangle.Gauss12;
    case 7
      tbl = triangle.Gauss16;
    case 8
      tbl = triangle.Gauss16;
    case 9
      tbl = triangle.Gauss19;
    case 10
      tbl = triangle.Gauss25;
    case 11
      tbl = triangle.Gauss33;
    case 12
      tbl = triangle.Gauss33;
    case 13
      tbl = triangle.Gauss37;
    case 14
      tbl = triangle.Gauss42;
    otherwise
      error('quadrature with order %d not implemented',r);
  end
end