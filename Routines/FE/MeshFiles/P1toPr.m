function [pr,tr] = P1toPr(p,t,r,barycentric)
% Author: Jet Hoe Tang

% Nodes
A = p(t(:,1),:);
B = p(t(:,2),:);
C = p(t(:,3),:);

% Complete Nodes
xx = barycentric*[A(:,1)';B(:,1)';C(:,1)'];
yy = barycentric*[A(:,2)';B(:,2)';C(:,2)'];

% Sort: [S;M;G]
IS = (1:3);
IM = (1: (3*(r-1))) + 3;
IG = (3*r+1):size(barycentric,1);
plist  = [...
          reshape(xx(IS,:),[],1),reshape(yy(IS,:),[],1);...
          reshape(xx(IM,:),[],1),reshape(yy(IM,:),[],1);...
          reshape(xx(IG,:),[],1),reshape(yy(IG,:),[],1);...
          ];
elindex = reshape(1:size(plist,1),size(barycentric,1),size(t,1));
maps = inversemapping([...
                       reshape(elindex(IS,:),[],1);
                       reshape(elindex(IM,:),[],1);...
                       reshape(elindex(IG,:),[],1);...
                       ]);
elist = maps(elindex)';

% Remove duplicates
precisiondecimal = floor(log(eps)/log(.1))-4;
[~,expmap,redmap] = unique(round(plist,precisiondecimal),'rows');
filterednodes = plist(expmap,:);
smap = unique(expmap);
filter = zeros(size(redmap,1),1);
filter(smap) = 1:size(smap,1);
pr = filterednodes(redmap(smap),:);
tr = filter(expmap(redmap(elist)));
if size(t,1) == 1; tr = tr(:)'; end
end