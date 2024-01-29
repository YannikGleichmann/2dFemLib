function fig = PlotFEFunction(f,saveFlag,view2,varargin)
% This function plots the Finite Element Function for P1, P2, P3, P1b, P2b,
% P3b functions
% 
% -------------------------------------------------------------------------
% Usage:
%       fig = PlotFEFunction(f)
%       fig = PlotFEFunction(f,varargin)
% 
% -------------------------------------------------------------------------
% Input:
% * f            - argument of Function class
% * varargin     - varargin as in trisurf;
%                  type help trisurf for detailed explanation
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
if isempty(varargin) || ~exist('varargin','var')
    varargin = {'facecolor','interp','edgecolor','none'};
end
if ~exist('f','var') || ~isobject(f)
    error('first input must be a Function class');
end
if exist('saveFlag', 'var') == 0 || isempty(saveFlag); saveFlag = 0; end
if exist('view2', 'var') == 0 || isempty(view2); view2 = 1; end
% deterime FunctionSpace for f and choose correct asis
[r,b] = FunctionSpace2Dim(f.Mesh);
% since trisurf can only handle P1 meshes, every triangle better than P1
% (P2, ...) has to be split up in subtriangles
if r == 1
    % P1 or P1b; basically here is nothing to do
    PlotMatrix = [1,2,3];
elseif b
    r = r-1;    % recall that for Pkb, the dimension is r = k+1 (not k!)
    if r == 2
        % P2b
        PlotMatrix = [1,4,7;
                      1,6,7;
                      3,5,6;
                      2,4,5;
                      4,5,7;
                      5,6,7];
    elseif r == 3
        % P3b
        PlotMatrix = [1,4,10;
                      1,9,10;
                      2,5,6;
                      3,7,8;
                      4,5,11;
                      4,10,11;
                      5,6,11;
                      6,7,11;
                      7,8,12;
                      7,11,12;
                      8,9,12;
                      9,10,12;
                      10,11,12];
    else
        error('not implemented');
    end
elseif ~b
    if r == 2
        % P2
        PlotMatrix = [1,4,6;
                      2,4,5;
                      3,5,6;
                      4,5,6];
    elseif r == 3
        % P3
        PlotMatrix = [1,9,10;
                      1,4,10;
                      2,5,6;
                      3,7,8;
                      4,5,6;
                      4,6,10;
                      6,7,10;
                      7,9,10;
                      7,8,9];
    else
        error('not implemented');
    end
else
    error('FunctionSpace not recognized');
end

p = f.Mesh.p;
t = f.Mesh.t;
t = reshape(t(:,PlotMatrix')',3,[])';
if any(size(f.NodalValues) == 1)
    fig = trisurf(t,p(:,1),p(:,2),f.NodalValues,...
         varargin{:});
else
    [~,MaxIt] = size(f.NodalValues);
    F = f.NodalValues;
    mi = min(F(:)); ma = max(F(:));
    if saveFlag
        disp('start movie; no visibility');
        set(0,'DefaultFigureVisible','off');
        writerObj = VideoWriter('VideoFile','MPEG-4');
        writerObj.FrameRate = 30;
        open(writerObj);
        for k = 1:MaxIt
            fig = trisurf(t,p(:,1),p(:,2),F(:,k),varargin{:});
            zlim([mi,ma]);
            if view2; view(2); axis equal; axis tight; end
            writeVideo(writerObj,getframe(gcf));
        end
        close(writerObj);
        set(0,'DefaultFigureVisible','on');
        close all;
        disp('done saving movie file');
    else
        for k = 1:MaxIt
        fig = trisurf(t,p(:,1),p(:,2),F(:,k),varargin{:});
        zlim([mi,ma]);
        if view2; view(2); axis equal; axis tight; end
        drawnow;
        end
    end
end