function PlotQuadratureNodes(tbl)
% this function plots the quadrature nodes obtained from tbl, see the
% function QuadratureCoefficientTriangle for more information

T = [0,0;1,0;0,1;0,0];
plot(T(:,1),T(:,2),'b',tbl(1,:),tbl(2,:),'xr');