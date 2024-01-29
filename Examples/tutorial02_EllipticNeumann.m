%% START
% This code solves the ellpitic problem with zero neumann conditions
%           - u'' + u = f on X := (0, 1)^2
%               du/dn = 0 on dX (boundary)
% with solution u given as u(x, y) = (2*pi*x - sin(2*pi*x)).*cos(2*pi*y)

%% PATH
clear variables; close all; clc;
restoredefaultpath
addpath(genpath('../Routines'));

%% exact solution, right hand sinde
u_   = @(x,y) (2*pi*x - sin(2*pi*x)).*cos(2*pi*y);
ux_  = @(x,y) 2*pi*(1 - cos(2*pi*x)).*cos(2*pi*y);
uy_  = @(x,y) (2*pi*x - sin(2*pi*x)).*sin(2*pi*y)*(-2*pi);
uxx_ = @(x,y) (2*pi)^2*sin(2*pi*x).*cos(2*pi*y);
uyy_ = @(x,y) (2*pi*x - sin(2*pi*x)).*cos(2*pi*y)*(-(2*pi)^2);
f_   = @(x,y) - uxx_(x,y) - uyy_(x,y) + u_(x,y);

%% function space, mesh
% for convinience use the same function space and mesh for all functions
FunctionSpace = 'P3';
mesh = genMesh('hmax',0.05,'FunctionSpace',FunctionSpace);

%% Functions
uExact = Function(mesh,u_(mesh.p(:,1),mesh.p(:,2)));
f = Function(mesh,f_(mesh.p(:,1),mesh.p(:,2)));

%% FE matricies
A = StiffnessMatrix2D(mesh);
M = MassMatrix2D(mesh);
L = LoadVector2D(mesh,f);

%% Compute FE solution
U = full((A+M)\L);
u = Function(mesh,U);

%% Plot
figure(1);
subplot(1,2,1); u.Plot; axis square; title('FE solution');
subplot(1,2,2); uExact.Plot; axis square; title('exact solution');