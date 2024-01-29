function tbl = Order2Quadrature1D(r,b)
% this function return a quadrature rule of order r as a table tbl where
% row 1 corresponds to the x value and row 2 are the weights to the 
% corresponding quadrature points.
% 
% -------------------------------------------------------------------------
% Usage:
%       tbl = Order2Quadrature1D(r)
%       tbl = Order2Quadrature(r,[]) ( = Order2Quadrature(r) )
%       tbl = Order2Quadrature1D(r,b)
% 
% -------------------------------------------------------------------------
% Input:
%   * r     - order of Gaussian quadrature rule
%   * b     - =1 if quadrature rule for mass lumping is desired
% 
% -------------------------------------------------------------------------
% Output:
%   * tbl = xi(1),   xi(2), ..., xi(nQuad)      <- x values
%           w(1),     w(2), ..., w(nQuad)       <- weights
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   14.06.2020

if exist('b','var') == 0 || isempty(b); b = 0; end
interval = QuadratureCoefficientInterval();

if b
  switch r
    case 1
      tbl = interval.Bubble2;
    case 2
      error('no bubble space with dimension %d',r');
    case 3
      tbl = interval.Bubble3;
    case 4
      tbl = interval.Bubble4;
    otherwise
      error('not implemented');
  end
else
  switch r
    case 1
      tbl = interval.Gauss1;
    case 2
      tbl = interval.Gauss1;
    case 3
      tbl = interval.Gauss2;
    case 4
      tbl = interval.Gauss2;
    case 5
      tbl = interval.Gauss3;
    case 6
      tbl = interval.Gauss3;
    case 7
      tbl = interval.Gauss4;
    case 8
      tbl = interval.Gauss4;
    case 9
      tbl = interval.Gauss5;
    case 10
      tbl = interval.Gauss5;
    case 11
      tbl = interval.Gauss6;
    case 12
      tbl = interval.Gauss6;
    case 13
      tbl = interval.Gauss7;
    case 14
      tbl = interval.Gauss7;
    otherwise
      error('quadrature with order %d not implemented',r);
  end
end