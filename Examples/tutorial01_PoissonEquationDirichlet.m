%% START
% This example solves the poisson equation with variation and zero diriclet
% boundary conditions. Consider
%           -(a u')' = f in X := (0,1)^2
%                  u = 0 on dX (boundary)
% Let the solution be given as u(x,y) = sin(k*pi*x) .* sin(k*pi*y) and the
% variation as  x^2 + y^2

%% PATH
clear variables; close all; clc;
restoredefaultpath;
addpath(genpath('../Routines'));

%% variation, exact solution, right hand side (rhs)
% Define the variation, exact solution, and right hand side f as function
% handles. These will be indicated with an underscore at the end
k    = 2;
u_   = @(x,y) sin(k*pi*x) .* sin(k*pi*y);
a_   = @(x,y) x.^2 + y.^2;
ux_  = @(x,y) k*pi*cos(k*pi*x) .* sin(k*pi*y);
uy_  = @(x,y) k*pi*sin(k*pi*x) .* cos(k*pi*y);
uxx_ = @(x,y) -k^2*pi^2*sin(k*pi*x) .* sin(k*pi*y);
uyy_ = @(x,y) -k^2*pi^2*sin(k*pi*x) .* sin(k*pi*y);
ax_  = @(x,y) 2*x;
ay_  = @(x,y) 2*y;
f_   = @(x, y) -(ax_(x,y).*ux_(x,y) + ay_(x,y).*uy_(x,y) ...
               + a_(x,y).*uxx_(x,y) + a_(x,y).*uyy_(x,y));
         
%% function spaces, Mesh (class)
% Next, define the function space for the solution (FunctionSpace), the
% variation (FunctionSpaceVar) and the rhs (FunctionSpaceRhs). Then
% generate the different meshes for those with genMesh.
FunctionSpace    = 'P3';
FunctionSpaceVar = 'P2';
FunctionSpaceRhs = 'P3';
mesh    = genMesh('hmax',0.05,'FunctionSpace',FunctionSpace);
meshVar = genMesh('hmax',0.05,'FunctionSpace',FunctionSpaceVar);
meshRhs = genMesh('hmax',0.05,'FunctionSpace',FunctionSpaceRhs);

%% Function (class)
% Now, define a Function class for every function used and evaluate the
% function handle on the given mesh for this function. Lower case variables
% are of Function class, upper case variables are for the nodal values of
% the function on the grid points
uExact = Function(mesh,u_(mesh.p(:,1),mesh.p(:,2)));
a = Function(meshVar,a_(meshVar.p(:,1),meshVar.p(:,2)));
f = Function(meshRhs,f_(meshRhs.p(:,1),meshRhs.p(:,2)));

%% FE matrix
% To compute the FE solution we need the stiffness matrix with variation a
% and the load vector. Here the variation a and the rhs f should be of
% Function class. It can also be a function handle (then it will be thread
% as it is a Function lying in mesh or it can be an (numerical) array (then
% it will be thread as the nodal values of it at mesh
A = StiffnessMatrix2D(mesh,a);
L = LoadVector2D(mesh,f);
% Now we have to incorporate the zero Dirichlet boundary conditions.
% Therefore we use the Mesh function idxBoundaryInteriorEdges which gives
% us the index of the boundary and interior points of the mesh class
[bnd,int] = idxBoundaryInteriorEdges(mesh);

%% Compute the solution
% Compute the FE solution in the interior (nodal values), incorporat the
% boundary conditions put it in the Function class
U = zeros(numberDof(mesh),1);
U(int) = A(int,int)\L(int);
u = Function(mesh,U);

%% Plot
% Plot the FE solution and the exact solution with the PlotFEFunction
% routine in a subplot and compare them with the eye-ball-metric
figure(1);
subplot(1,2,1); u.Plot; axis square; title('FE solution');
subplot(1,2,2); uExact.Plot; axis square; title('exact solution');