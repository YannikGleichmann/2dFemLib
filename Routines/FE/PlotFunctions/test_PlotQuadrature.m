clear; close all; clc;
% add all (sub)folder in current m-Path to search Path %
folder = fileparts(which(mfilename));                  %
addpath(genpath(folder));                              %
%------------------------------------------------------%
% triangle = QuadratureCoefficientTriangle();
% PlotQuadratureNodes(triangle.Gauss27);

tbl = Order2Quadrature(2,1);
PlotQuadratureNodes(tbl);