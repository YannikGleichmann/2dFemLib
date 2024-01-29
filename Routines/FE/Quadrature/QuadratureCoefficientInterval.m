function interval = QuadratureCoefficientInterval()

% The function return a struct with the different quadrature formulas for
% the 2D unit triangle T = (0,0), (1,0), (0,1).
%
% The quadrature formulas are represented by a 3 x nQuadNodes array
% * row 1: xi       represents the cartesian coordinate x of a quadrature
%                   point
% * row 3:  w       stores the weight of the corresponding quadrature
%                   point. The weights are scaled by
%                   sum(weights) = volume of the unit interval = 1;
% 
% Reference:
% * https://de.wikipedia.org/wiki/Gauß-Quadratur
% * https://pomax.github.io/bezierinfo/legendre-gauss.html
% 
% Author: Yannik Gleichmann
% Date: 22.06.2020

t = @(s) (s+1)/2;
%% Gauss-Legendre
% Gauss 2nd order
tbl = [0,2];
Gauss1 = [t(tbl(:,1))'; tbl(:,2)'/sum(tbl(:,2))];
% Gauss 4th order
tbl = [-sqrt(1/3),1;
        sqrt(1/3),1];
Gauss2 = [t(tbl(:,1))'; tbl(:,2)'/sum(tbl(:,2))];
% Gauss 6th order
tbl = [-sqrt(3/5), 5/9;
                0, 8/9;
        sqrt(3/5), 5/9];
Gauss3 = [t(tbl(:,1))'; tbl(:,2)'/sum(tbl(:,2))];
% Gauss 8th order
tbl = [-sqrt(3/7+2/7*sqrt(6/5)), (18-sqrt(30))/36;
       -sqrt(3/7-2/7*sqrt(6/5)), (18+sqrt(30))/36;
        sqrt(3/7-2/7*sqrt(6/5)), (18+sqrt(30))/36;
       -sqrt(3/7+2/7*sqrt(6/5)), (18-sqrt(30))/36];
Gauss4 = [t(tbl(:,1))'; tbl(:,2)'/sum(tbl(:,2))];
% Gauss 10th order
tbl = [-0.9061798459386640, 0.2369268850561891;
       -0.5384693101056831, 0.4786286704993665;
        0.0000000000000000, 0.5688888888888889;
        0.5384693101056831, 0.4786286704993665;
        0.9061798459386640, 0.2369268850561891];
Gauss5 = [t(tbl(:,1))'; tbl(:,2)'/sum(tbl(:,2))];
% Gauss 12th order
tbl = [-0.9324695142031521, 0.1713244923791704;
       -0.6612093864662645, 0.3607615730481386;
       -0.2386191860831969, 0.4679139345726910;
        0.2386191860831969, 0.4679139345726910;
        0.6612093864662645, 0.3607615730481386;
        0.9324695142031521, 0.1713244923791704];
Gauss6 = [t(tbl(:,1))'; tbl(:,2)'/sum(tbl(:,2))];
% Gauss 14th order
tbl = [-0.9491079123427585, 0.1294849661688697;
       -0.7415311855993945, 0.2797053914892766;
       -0.4058451513773972, 0.3818300505051189;
        0.0000000000000000, 0.4179591836734694;
        0.4058451513773972, 0.3818300505051189;
        0.7415311855993945, 0.2797053914892766;
        0.9491079123427585, 0.1294849661688697];
Gauss7 = [t(tbl(:,1))'; tbl(:,2)'/sum(tbl(:,2))];

%% mass lumping quadrature
% P1b
tbl = [0,1;
       1,1];
Bubble2 = [tbl(:,1)'; tbl(:,2)'/sum(tbl(:,2))];
% P2b
tbl = [0, 1/6;
       0.5,4/6;
       1,1/6];
Bubble3 = [tbl(:,1)'; tbl(:,2)'/sum(tbl(:,2))];
% P3b
tbl = [0,0.293469555909040,0.706530444090960,1;
       0.098093695372309,0.401906304627690,0.401906304627692,0.098093695372308];
Bubble4 = [tbl(1,:); tbl(2,:)/sum(tbl(2,:))];
%% Summary
interval = struct();

% Gauss
interval.Gauss1 = Gauss1;
interval.Gauss2 = Gauss2;
interval.Gauss3 = Gauss3;
interval.Gauss4 = Gauss4;
interval.Gauss5 = Gauss5;
interval.Gauss6 = Gauss6;
interval.Gauss7 = Gauss7;

% Bubble
interval.Bubble2 = Bubble2;
interval.Bubble3 = Bubble3;
interval.Bubble4 = Bubble4;
