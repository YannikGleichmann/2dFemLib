function ID = PointInTriangle(P,t,p)

if ~(size(P) == 2); error('Wrong point matrix dimension'); end
P = P';

A = p(t(:,1),:)';
B = p(t(:,2),:)';
C = p(t(:,3),:)';

nP = size(P,2);
nTri = size(t,1);

iTri = mod(0:nP*nTri-1,nTri)+1;
iP = ((1:nP*nTri) - iTri)/nTri + 1;

v0 = C(:,iTri) - A(:,iTri);
v1 = B(:,iTri) - A(:,iTri);
v2 = P(:,iP) - A(:,iTri);

% Compute dot products
dot00 = dot(v0,v0);
dot01 = dot(v0,v1);
dot02 = dot(v0,v2);
dot11 = dot(v1,v1);
dot12 = dot(v1,v2);

% Compute barycentric coordinates
invArea = 1./(dot00.*dot11-dot01.*dot01);
u = (dot11.*dot02-dot01.*dot12).*invArea;
v = (dot00.*dot12 - dot01.*dot02).*invArea;

boo = sparse(1,nTri*nP);
boo( (u >= 0) & (v >= 0) & (u + v <= 1)) = 1;
boo = reshape(boo,nTri,nP)';
[~,ID] = max(boo,[],2);