clear; close all; clc;
x_min = 0; x_max = 1; y_min = 0; y_max = 1;
r = RectangleGeometry(x_min, x_max, y_min, y_max);
[p, e, t] = initmesh(r, 'hmax', 0.6, 'Jiggle', 'mean');
mesh = Mesh;
mesh.g = r;     % geometry points
mesh.p = p;     % point list
mesh.e = e;     % edge list
mesh.t = t;     % element/triangle list
figure(1);
pdemesh(mesh.p, mesh.e, mesh.t, 'NodeLabels', 'on', 'ElementLabels', 'on');

p = mesh.p;                 % point list
t = mesh.t;                 % element/triangle list
e = mesh.e;                 % triangle list

nDof = numberDof(mesh);     % number of points
nTri = numberTri(mesh);     % number of elements/triangles

p = p';
t = t(1:3, :)';
f = @(x, y) x + y;
fNodes = f(p(:, 1), p(:, 2));

p1 = rand(5, 2);
load('random_points.mat')
% p1 = p;
fInter = zeros(size(p1, 1), 1);

idx_binned_points = ones(size(p1, 1), 1);
p1_copy = p1;
N = {@(x, y) 1 - x - y, @(x, y) x, @(x, y) y};
nBasis = length(N);
for k = 1:nTri
    for i = 1:size(p1_copy, 1)
        poly = [p(t(k,1), :); p(t(k,2), :); p(t(k,3), :)];
        p_glob = p1_copy(i, :)';
        in = inpolygon(p_glob(1), p_glob(2), poly(:, 1), poly(:, 2));
        if in
            p1_copy(i, :) = [inf, inf];
            P1 = p(t(k,1),:)';
            P2 = p(t(k,2),:)';
            P3 = p(t(k,3),:)';
            J_K = [P2 - P1, P3 - P1];
            p_loc = J_K\p_glob - J_K\P1;
            for b = 1:nBasis
                fInter(i) = fInter(i) + fNodes(t(k, b))*N{b}(p_loc(1), p_loc(2));
            end
        end
    end
end
fInter == f(p1(:, 1), p1(:, 2));