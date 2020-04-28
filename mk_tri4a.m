function [tri,xyz,nt,slx] = mk_tri4a(dat,tol,iplt);
%MK_TRI4A Makes a triangular mesh by using the ordered slice data from
%         the digitized MRI and using angles from polar coordinates.
%
%         [TRI,XYZ,NT] = MK_TRI4A(DAT) given a cell array containing
%         three (3) columns matrices with slice coordinate point data,
%         DAT, returns the three (3) column triangle connectivity
%         matrix, TRI.  The coordinates for the patella are returned in
%         a three (3) column matrix, XYZ.  The number of returned
%         triangles, NT, may also be returned.
%
%         NOTES:  1.  Each slice coordinate data matrix must correspond
%                 to one index into the cell array DAT.
%
%                 2.  The angles along each slice are used to
%                 determine the triangulation.
%
%                 3.  The M-files plane_fit.m, sl_info.m, tri_fix2.m and
%                 tri_norm.m must be in the current path or directory.
%
%         12-Oct-2015 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  iplt = false;
end
%
if iplt
  hf = figure;
  orient tall;
end
%
if (nargin<2)
  tol = 0.1;            % Tolerance on dot product (projection of directions)
end
%
if (nargin<1)
  error(' *** ERROR in MK_TRI4A:  No input data!');
end
%
% Get Number of Slices
%
[nslice,npts,n] = sl_info(dat);
tri = [];
slx = zeros(nslice-1,1);
%
% Loop through Pairs of Slices
%
for k = 2:nslice
%
% Get Slices into a Plane
%
   xyz1 = dat{k-1};     % Coordinates for first slice
   xyz2 = dat{k};       % Coordinates for second slice
   npt1 = npts(k-1);    % Number of points in first slice
   npt2 = npts(k);      % Number of points in second slice
   cntr = mean([xyz1; xyz2]);          % Center of slices
   [~,~,~,r] = plane_fit(xyz1(:,1),xyz1(:,2),xyz1(:,3));   % Rotation matrix
%
% Transform Slice Data
%
   xyz1t = (xyz1-repmat(cntr,npt1,1))*r;
   xyz2t = (xyz2-repmat(cntr,npt2,1))*r;
   slx(k-1) = mean(xyz2t(:,3))-mean(xyz1t(:,3));
%
% Get Polar Coordinates (Angles in radians)
%
   [ang1,rad1] = cart2pol(xyz1t(:,1),xyz1t(:,2));
   idx = find(ang1>pi/2);              % Check for jump at +pi/2
   ang1(idx) = ang1(idx)-2*pi;
%
   [ang2,rad2] = cart2pol(xyz2t(:,1),xyz2t(:,2));
   idx = find(ang2>pi/2);              % Check for jump at +pi/2
   ang2(idx) = ang2(idx)-2*pi;
%
   ang2 = ang2-ang1(1);
   ang1 = ang1-ang1(1);
%
% Get Directions of Angles
%
   dang = diff(ang1);
   dang = sum(dang);
   dir1 = sign(dang);
%
   dang = diff(ang2);
   dang = sum(dang);
   dir2 = sign(dang);
%
% Make Sure Angles are in Order
%
   if dir1>0
     ang1 = sort(ang1,1,'ascend');
   else
     ang1 = sort(ang1,1,'descend');
   end
%
   if dir2>0
     ang2 = sort(ang2,1,'ascend');
   else
     ang2 = sort(ang2,1,'descend');
   end
%
% Delaunay Triangulation
%
   xt = [zeros(npt1,1); slx(k-1)*ones(npt2,1)];
   yt = [ang1; ang2];
   tril = delaunay(xt,yt);
%
   if iplt
     ntril = size(tril,1);
     cla;
     plot(xt,yt,'k.');
     hold on;
     trimesh(tril,xt,yt);
     text(xt,yt,int2str((1:length(xt))'),'Color','b','FontSize',12);
     text(mean(xt(tril),2),mean(yt(tril),2),int2str((1:ntril)'), ...
          'Color','r','FontSize',12);
     pause;
   end
%
% Improve Mesh
%
   xyz = [xyz1; xyz2];
   tril = tri_fix2(tril,xyz);
%
% Check Normals
%
   [nx,ny,nz,xc,yc,zc] = tri_norm(tril,xyz);
   ntril = size(tril,1);
   nv = [nx ny nz];
   xyzc = [xc yc zc]-repmat(cntr,ntril,1);
   irev = false(ntril,1);
   for l = 1:ntril
      irev(l) = xyzc(l,:)*nv(l,:)'<0;
   end
   if nnz(irev)>ntril/2 % If most are reversed, then all should be reversed
     tril = tril(:,[1 3 2]);
   end
%
% Global Node IDs
%
   nid = n(k-1)+1:n(k+1);
   tri = [tri; nid(tril)];
%
end
%
% Improve Mesh
%
xyz = cell2mat(dat);
tri = tri_fix2(tri,xyz);
%
nt = size(tri,1);
%
return