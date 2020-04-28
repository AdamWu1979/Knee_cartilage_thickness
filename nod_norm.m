function [nodv,arean,triv,areat,cg] = nod_norm(tri,xyz,iplt);
%NOD_NORM  Calculates corner point (node) normals of a triangular
%          surface mesh. 
%
%          NODV = NOD_NORM(TRI,XYZ) returns a matrix of normal vectors
%          at the corner points (nodes) of a triangular surface mesh in
%          the three (3) column matrix, NODV.  The triangular surface
%          mesh is defined by a three (3) column triangle connectivity
%          matrix, TRI, and the X, Y and Z coordinates of the nodes in
%          a three (3) column matrix, XYZ.
%
%          [NODV,AREAN,TRIV,AREAT,CG] = NOD_NORM(TRI,XYZ) returns the
%          areas of the triangles connected to the node in the vector,
%          AREAN, normal vectors of the triangles in a three (3) column
%          matrix, TRIV, the areas of the triangles in the vector,
%          AREAT, and the X, Y and Z coordinates of the area centers of
%          the triangles in the three (3) column matrix CG.
%
%          NOTES:  None.
%
%          18-April-2011 * Mack Gardner-Morse
%

%#######################################################################
%
% Check Input Parameters
%
if (nargin<2)
  error(' *** ERROR in nod_norm:  Not enough input arguments.');
end
%
if (nargin<3)
  iplt = false;
end
%
% Check Inputs
%
ncol1 = size(tri,2);
ncol2 = size(xyz,2);
%
if (ncol1~=3)&(ncol2~=3)
  error([' *** ERROR in nod_norm:  Input matrices must have three', ...
        ' (3) columns!']);
end
%
% Find Triangle Normals and Areas
%
nt = size(tri,1);       % Number of triangles
nnod = size(xyz,1);     % Number of nodes
%
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
%
xe = x(tri);
ye = y(tri);
ze = z(tri);
%
if nt==1
  xe = xe';
  ye = ye';
  ze = ze';
end
%
v1 = [xe(:,2)-xe(:,1) ye(:,2)-ye(:,1) ze(:,2)-ze(:,1)]';
v2 = [xe(:,3)-xe(:,1) ye(:,3)-ye(:,1) ze(:,3)-ze(:,1)]';
%
triv = cross(v1,v2);
%
areat = sqrt(sum(triv.*triv))';        % Norm of normal vectors
rnorm = repmat(areat,1,3);
%
triv = triv'./rnorm;    % Normalize normal vectors
%
areat = areat/2;        % Area of triangles
rnorm = repmat(areat,1,3);             % Duplicate so area can be used as a weight
%
% Find Node Normals
%
nodv = zeros(nnod,3);
arean = zeros(nnod,1);
%
atriv = rnorm.*triv;    % Weight by area of triangles
%
for k = 1:nt
   idn = tri(k,:)';
   nodv(idn,:) = nodv(idn,:)+repmat(atriv(k,:),3,1);
   arean(idn) = arean(idn)+repmat(areat(k),3,1);
end
%
idg = find(arean);      % Find nodes with areas to prevent division by zero
nodv(idg,:) = nodv(idg,:)./repmat(arean(idg),1,3);    % Divide by total area (not necessary?)
%
rnorm = sqrt(sum(nodv(idg,:).^2,2));
nodv(idg,:) = nodv(idg,:)./repmat(rnorm,1,3);
%
% Find Centers
%
if nargout>4
  cg = [mean(xe,2) mean(ye,2) mean(ze,2)];
end
%
return