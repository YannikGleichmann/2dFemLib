%% Start
% Solves the wave equation with zero dirichlet boundary conditions
%       y_tt - y'' = f   in X x I := (0, 1)^2 x [0, T]
%          y(-, 0) = g_1 in X
%        y_t(-, 0) = g_2 in X
%                y = g_3 on dX x I
% where the solution is given as
% y(x, y, t) = sin(c*pi*x)*sin(c*pi*y)*sin(4*pi*t)

%% PATH
clear variables; close all; clc;
restoredefaultpath
addpath(genpath('../Routines'));

%% exact solution
c = 2; a = 2;
y_   = @(x,y,t) sin(c*pi*x).*sin(c*pi*y).*sin(a*pi*t);
yt_  = @(x,y,t) a*pi*sin(c*pi*x).*sin(c*pi*y).*cos(a*pi*t);
ytt_ = @(x,y,t) -(a*pi).^2*y_(x,y,t);
yxx_ = @(x,y,t) -(c*pi).^2*y_(x,y,t);
yyy_ = @(x,y,t) -(c*pi).^2*y_(x,y,t);

%% given data
f_  = @(x,y,t) ytt_(x,y,t) - (yxx_(x,y,t)+yyy_(x,y,t));
g1_ = @(x,y) y_(x,y,0);
g2_ = @(x,y) yt_(x,y,0);
g3_ = @(x,y,t) y_(x,y,t);

%% function space, mesh
FunctionSpace = 'P3b';
mesh = genMesh('hmax',0.05,'FunctionSpace',FunctionSpace);

%% FE maticies
% For mass lumping (indicated by P1b, P2b, P3b for the function space), the
% mass matrix is diagonal and can be inverted easily. To make the code more
% general we also allow the "normal" elements.
% MassMatrix2D will automatically diagonalize the mass matrix, if the mesh
% contains bubble functions. To surpess this, use
% MassMatrix2D(mesh,m,0) or
% MassMatrix2D(mesh,[],0) if no variation m is needed. If the last
% parameter is set to 0, the mass matrix will not be lumped. If it is set
% to 1 (or omited), the mass matrix will be lumped if an apropriate
% FunctionSpace is used. For more information see description of
% MassMatrix2D().
M = MassMatrix2D(mesh);
isdiagM = isdiag(M);
A = StiffnessMatrix2D(mesh);
if isdiag(M)
  Minv = spdiags(1./spdiags(M),0,numberDof(mesh),numberDof(mesh));
end

%% time
% Since we are working in time, we have to use the time class which stores
% the initial and end time, the time step length dt, and the number of time
% steps to reach from the initial to the end time. 
time    = Time;
time.T0 = 0;
time.T1 = 1;
m       = eigs(M,1,'sm','Tolerance',1e-3,'MaxIterations',350);
a       = eigs(round((A+A')/2,10),1,'lm','Tolerance',1e-3,...
               'MaxIterations',350);
dt      = 0.8*sqrt(2*m/a);
nT      = ceil((time.T1 - time.T0)/dt);
dt      = (time.T1 - time.T0)/nT;
time.nT = nT;
time.dt = dt;
t = linspace(time.T0,time.T1,time.nT+1);

%% source, boundary condition
F  = f_(mesh.p(:,1),mesh.p(:,2),t);
G1 = g1_(mesh.p(:,1),mesh.p(:,2));
G2 = g2_(mesh.p(:,1),mesh.p(:,2));
G3 = g3_(mesh.p(:,1),mesh.p(:,2),t);
f  = Function(mesh,F);
g1 = Function(mesh,G1);
g2 = Function(mesh,G2);
g3 = Function(mesh,G3);

%% initial conditions for leap frog
[bnd,int] = idxBoundaryInteriorEdges(mesh);
M = M(int,int);
A = A(int,int);
Minv = Minv(int,int);
Y = zeros(length(int),nT+1);
y0 = G1(int);
if isdiagM
  y1 = Minv*(M*G1(int) + dt*M*G2(int) + dt^2/2*M*F(int,1) ...
             - dt^2/2*A*G1(int));
else
  y1 = M\(M*G1(int) + dt*M*G2(int) + dt^2/2*M*F(int,1)...
          - dt^2/2*A*G1(int));
end
Y(:,1) = y0;
Y(:,2) = y1;

%% Compute the FE solution
if isdiagM
  for k = 1:nT-1
    y2 = Minv*(2*M*y1 - M*y0 - dt^2*A*y1 + dt^2*M*F(int,k+1));
    y0 = y1;
    y1 = y2;
    Y(:,k+2) = y2;
  end
else
  for k = 1:nT-1
    y2 = M\(2*M*y1 - M*y0 - dt^2*A*y1 + dt^2*M*F(int,k+1));
    y0 = y1;
    y1 = y2;
    Y(:,k+2) = y2;
  end
end

Ytemp = G3; Ytemp(int,:) = Y;
%% Plot
ySnip = Function(mesh,Ytemp(:,round(nT/4)));
yExactSnip = Function(mesh,G3(:,round(nT/4)));
figure();
subplot(1,2,1); ySnip.Plot; axis square;
        title(sprintf('FE solution at time t = %f sec',dt*round(nT/4)));
subplot(1,2,2); yExactSnip.Plot; axis square;
        title(sprintf('exact solution at time t = %f sec',dt*round(nT/4)));
% The f.Plot can also handle time dependent functions. For that the
% nodal values should be an numberDofs x nTimeSteps array. Note that this
% only plots the P1 equivalent.
y = Function(mesh,G3);
y.NodalValues(int,:) = Y;
pause(0.01);
y.Plot