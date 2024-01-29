function fig = PlotFEMesh(mesh,PlotNodes,PlotSubmesh)
% This function plots the Finite Element mesh for P1, P2, P3, P1b, P2b, P3b
% meshes
% 
% -------------------------------------------------------------------------
% Usage:
%       fig = PlotFEMesh(mesh)
%       fig = PlotFEMesh(mesh,PlotSubmesh)
% 
% -------------------------------------------------------------------------
% Input:
% * mesh            - argument of Mesh class
% * PlotSubmesh     - PlotFEMesh only plots the "P1" mesh;
%                     if PLOTSubmesh = 1:
%                     the outer triangle will be split up into subtriangles
%                     defines by the point matrix mesh.p and plotted in red
%                     dotted lines;
% 
% -------------------------------------------------------------------------
% Output:
% * patch object
% 
%--------------------------------------------------------------------------
% Author: Yannik Gleichmann
% Date:   17.06.2020

% folder = fileparts(which(mfilename));
% addpath(genpath(folder));
if ~exist('mesh','var') || ~isobject(mesh)
    error('first input must be of Mesh class');
end
if ~exist('PlotSubmesh','var') || isempty('PlotSubmesh')
    PlotSubmesh = 0;
end
if ~exist('PlotNodes','var') || isempty('PlotNodes')
    PlotNodes = 0;
end
% deterime FunctionSpace for f and choose correct asis
[r,b] = FunctionSpace2Dim(mesh);
% since trisurf can only handle P1 meshes, every triangle better than P1
% (P2, ...) has to be split up in subtriangles
if r == 1
    % P1 or P1b; basically here is nothing to do
    PlotMatrix = [1,2,3];
elseif b
    r = r-1;    % recall that for Pkb, the dimension is r = k+1 (not k!)
    if r == 2
        % P2b
        PlotMatrix = [1,4,7; 1,6,7; 3,5,6; 2,4,5; 4,5,7; 5,6,7];
    elseif r == 3
        % P3b
        PlotMatrix = [1,4,10; 1,9,10; 2,5,6; 3,7,8; 4,5,11; 4,10,11;
                      5,6,11; 6,7,11; 7,8,12; 7,11,12; 8,9,12;
                      9,10,12; 10,11,12];
    else
        error('not implemented');
    end
elseif ~b
    if r == 2
        % P2
        PlotMatrix = [1,4,6; 2,4,5; 3,5,6; 4,5,6];
    elseif r == 3
        % P3
        PlotMatrix = [1,9,10; 1,4,10; 2,5,6; 3,7,8; 4,5,6; 4,6,10;
                      6,7,10; 7,9,10; 7,8,9];
    else
        error('not implemented');
    end
else
    error('FunctionSpace not recognized');
end
p = mesh.p;
t = mesh.t;
fig = patch('Faces',t(:,1:3),'Vertices',p,'FaceColor','none',...
            'EdgeColor','black','LineWidth',1); hold on;
t = reshape(t(:,PlotMatrix')',3,[])';
P1 = p(t(:,1),:); P2 = p(t(:,2),:); P3 = p(t(:,3),:);
v1 = P2-P1; v2 = P3-P1;
idx = (v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2))<0;
t(idx,:) = t(idx,[2,1,3]);
if PlotSubmesh && (r~=1)
    patch('Faces',t, ...
           'Vertices',p,'FaceColor','none','EdgeColor','red',...
           'Linestyle',':');
end
if PlotNodes
    plot(mesh.p(:,1),mesh.p(:,2),'xr'); hold off;
end
hold off; axis equal;

% pdemesh(mesh.p',mesh.t','NodeLabels','on','ElementLabels','on')