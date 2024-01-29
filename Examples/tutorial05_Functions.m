%% START
% This example will explain the Function class and its properties in detail

%% PATH
% The usual procedure; and add the path with the routines
clear; close all; clc; restoredefaultpath;
addpath(genpath('../Routines'));

%% MESH
FunctionSpace = 'P1';
mesh = genMesh('hmax',0.025,'FunctionSpace',FunctionSpace);

%% FUNCTION
f_ = @(x,y) sin(2*pi*x).*sin(2*pi*y);
F = f_(mesh.p(:,1),mesh.p(:,2));
f = Function(mesh,F);
figure();
subplot(3,2,1); f.Plot; title('f'); axis square;
% We can also perfom arithmetic operations such as +, -, *, /, ^n, sqrt
% to a function by using oop:
g = Function(mesh,cos(2*pi*mesh.p(:,1)).*cos(2*pi*mesh.p(:,2)));
h_plus = plus(f,g);
h_minus = minus(f,g);
h_times = times(f,g);
f_squared = power(f,2);
% and plot them
subplot(3,2,2); g.Plot; title('g'); axis square;
subplot(3,2,3); h_plus.Plot; title('f+g'); axis square;
subplot(3,2,4); h_minus.Plot; title('f-g'); axis square;
subplot(3,2,5); h_times.Plot; title('f*h'); axis square;
subplot(3,2,6); f_squared.Plot; title('f^2'); axis square;
% The next thing we can do is calculating varius norms and errors:
L2norm_f = L2norm(f);
semiH1norm_f = semiH1norm(f);
error = L2error(h_plus,h_minus);
fprintf(['L2-norm of f: %f\nsemi-H1-norm of f: %f\nL2-error ' ...
         'of h_plus and h_minus: %f\n'],L2norm_f,semiH1norm_f,error);