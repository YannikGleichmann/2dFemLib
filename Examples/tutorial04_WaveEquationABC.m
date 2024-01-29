%% Start
% Solves the wave equation with absorbing boundary conditions
%         y_tt - (u y')' = f in (0, 1)^2 x [0, T] =: O x I
%                 y(t=0) = g1 in O
%               y_t(t=0) = g2 in O
%            y_t + dy/dn = g3 on dO x I
% with leap frog, where the given solution is
%   y(x, y, t) = (2*pi*x - sin(2*pi*x))*cos(2*pi*y)*sin(t)
% hence dy/dn = 0 on dO and k = 0;
% the medium is chosen as u(x, y) = 1 + x + y

%% Path
clear; close all; clc;
restoredefaultpath;
addpath(genpath('../Routines'));

%% given data / exact solution
yExact_ = @(x, y, t) (2*pi*x - sin(2*pi*x)).*cos(2*pi*y).*sin(t);
u_      = @(x, y)     1 + x + y;
y_t_    = @(x, y, t) (2*pi*x - sin(2*pi*x)).*cos(2*pi*y).*cos(t);
y_tt_   = @(x, y, t)  -yExact_(x, y, t);
y_xx_   = @(x, y, t)  2*pi*cos(2*pi*y).*(1 - cos(2*pi*x) ...
                      + 2*pi*(1 + x + y).*sin(2*pi*x)).*sin(t);
y_yy_   = @(x, y, t)  2*pi*(2*pi*x - sin(2*pi*x)) ...
                      .*(-2*pi*(1 + x + y).*cos(2*pi*y) ...
                      - sin(2*pi*y)).*sin(t);
f_      = @(x, y, t)  y_tt_(x, y, t) - y_xx_(x, y, t) - y_yy_(x, y, t);

%% function space / mesh
FunctionSpace = 'P2b';
mesh = genMesh('hmax',0.05,'FunctionSpace',FunctionSpace);
FunctionSpaceVar = 'P1';
meshVar = genMesh('hmax',0.05,'FunctionSpace',FunctionSpaceVar);

%% Variation
% We need the variation first to determine the cfl condition for
% leap-frog-2
u = Function(meshVar,u_(meshVar.p(:,1),meshVar.p(:,2)));
usqrt = sqrt(u);
%% FE maticies
M    = MassMatrix2D(mesh);
Mbnd = MassBoundaryMatrix2D(mesh,usqrt);
A    = StiffnessMatrix2D(mesh,u);
isdiagM = isdiag(M) & isdiag(Mbnd);

if isdiagM
  nDof = mesh.numberDof;
  Minv = spdiags(1./spdiags(M),0,nDof,nDof);
end

%% time
time    = Time;
time.T0 = 0;
time.T1 = pi;
m       = eigs(M,1,'sm','Tolerance',1e-3,'MaxIterations',350);
a       = eigs(round((A+A')/2,10),1,'lm','Tolerance',1e-3,...
               'MaxIterations',350);
dt      = 0.9*sqrt(2*m/a);
nT      = ceil((time.T1 - time.T0)/dt);
dt      = (time.T1 - time.T0)/nT;
time.nT = nT;
time.dt = dt;
t = linspace(time.T0,time.T1,time.nT+1);

%% functions
yExact = Function(mesh,yExact_(mesh.p(:,1),mesh.p(:,2),t));
f = Function(mesh,f_(mesh.p(:,1),mesh.p(:,2),t));
F = f.NodalValues;
G1 = yExact_(mesh.p(:,1),mesh.p(:,2),0);
G2 = y_t_(mesh.p(:,1),mesh.p(:,2),0);
G3 = y_t_(mesh.p(:,1),mesh.p(:,2),t);
g1 = Function(mesh,G1);
g2 = Function(mesh,G2);
g3 = Function(mesh,G3);

%% initial conditions for leap frog
Y = zeros(numberDof(mesh),nT+1);
Y0 = G1;
if isdiagM
  Y1 = Minv*(M*G1 + dt*M*G2 + dt^2*M*F(:,1)/2 - dt^2*A*G1/2 ...
             + dt^2*Mbnd*G3(:,1)/2 - dt^2*Mbnd*G2/2);
else
  Y1 = M\(M*G1 + dt*M*G2 + dt^2*M*F(:,1)/2 - dt^2*A*G1/2 ...
          + dt^2*Mbnd*G3(:,1)/2 - dt^2*Mbnd*G2/2);
end
Y(:,1) = Y0;
Y(:,2) = Y1;

Ml = M+dt*Mbnd/2;
Mr = M-dt*Mbnd/2;
MA = 2*M - dt^2*A;

if isdiagM
  nDof = mesh.numberDof;
  Mlinv = spdiags(1./spdiags(Ml),0,nDof,nDof);
  for k = 1:nT-1
    Y2 = Mlinv*(MA*Y1 - Mr*Y0 + dt^2*M*F(:,k+1) + dt^2*Mbnd*G3(:,k+1));
    Y0 = Y1;
    Y1 = Y2;
    Y(:,k+2) = Y2;
  end
else
  for k = 1:nT-1
    Y2 = Ml\(MA*Y1 - Mr*Y0 + dt^2*M*F(:,k+1) + dt^2*Mbnd*G3(:,k+1));
    Y0 = Y1;
    Y1 = Y2;
    Y(:,k+2) = Y2;
  end
end

y = Function(mesh,Y);
y.Plot;