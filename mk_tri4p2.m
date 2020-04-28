function [tri,xyz,xyze,nt,slx] = mk_tri4p2(dat,iplt);
%MK_TRI4P2 Makes a triangular mesh by using the ordered slice data from
%         the digitized MRI patella data.
%
%         [TRI,XYZ,XYZE,NT] = MK_TRI4P2(DAT) given a cell array
%         containing three (3) columns matrices with slice coordinate
%         point data, DAT, returns the three (3) column triangle
%         connectivity matrix, TRI.  The coordinates for the patella
%         including two additional end points are returned in a three
%         (3) column matrix, XYZ.  The two additional end points are
%         returned in a two (2) row by three (3) column matrix, XYZE,
%         which are required to enclose the ends.  The number of
%         returned triangles, NT, may also be returned.
%
%         NOTES:  1.  Each slice coordinate data matrix must correspond
%                 to one index into the cell array DAT.
%
%                 2.  The coordinates should be ordered in the same
%                 direction in every slice.  The differences in angles
%                 in each slice are used to check the ordering
%                 direction.
%
%                 3.  The polar angle along each slice is used to
%                 determine the triangulation.  See mk_tri4.m and
%                 mk_tri4f.m for similar triangulations for the tibia
%                 and femur using arclengths.  See mk_tris.m for a
%                 triangulation based on each slice starting at the
%                 same in slice point.  See mk_triq.m for a quicker,
%                 but different triangulation based on the number of
%                 points in each slice.
%
%                 4.  The M-files in_tri2d.m, isect.m, plane_fit.m,
%                 sl_info.m, tri_fix2.m and tri_norm.m must be in the
%                 current path or directory.
%
%         24-Aug-2015 * Mack Gardner-Morse
%
%         23-Sep-2015 * Mack Gardner-Morse * Used minimum to find
%                                            matching angles.
%
%         30-Oct-2015 * Mack Gardner-Morse * Checks origin is within
%                                            both slices.
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  iplt = false;
end
%
if iplt
  hf = figure;
  orient tall;
end
%
if (nargin<1)
  error(' *** ERROR in MK_TRI4P2:  No input data!');
end
%
% Loop through Slices
%
dat = dat(:);
[nslice,npts] = sl_info(dat);
n = [0; cumsum(npts)];
%
tri = [];
slx = zeros(nslice-1,1);
%
for k = 2:nslice
%
% Get Center Point
%
   xyz1 = dat{k-1};
   xyz2 = dat{k};
%
   npt1 = npts(k-1);
   npt2 = npts(k);
%
   cntr = zeros(2,3);
   [cntr(1,:),nv,~,rot] = plane_fit(xyz1(:,1),xyz1(:,2),xyz1(:,3));
   cntr(2,:) = mean(xyz2);
   ctr = mean(cntr);    % Common center point
   cntr(1,:) = isect(cntr(1,:),nv,ctr,nv);       % Common center point within slice
   cntr(2,:) = isect(cntr(2,:),nv,ctr,nv);       % Common center point within slice
%
% Translate and Rotate Data
%
   xyz1t = xyz1-repmat(cntr(1,:),npts(k-1),1);
   xyz2t = xyz2-repmat(cntr(2,:),npts(k),1);
   xyz1t = xyz1t*rot;
   xyz2t = xyz2t*rot;
%
% Check that the Origin is Within the Slices
%
   tri1 = delaunay(xyz1t(:,1),xyz1t(:,2));
   in1 = in_tri2d(tri1,[xyz1t(:,1),xyz1t(:,2)],[0 0]);
   if ~in1
     ctr = mean(xyz1t);
     xyz1t = xyz1t-repmat(ctr,npts(k-1),1);
     xyz2t = xyz2t-repmat(ctr,npts(k),1);
   end
%
   tri2 = delaunay(xyz2t(:,1),xyz2t(:,2));
   in2 = in_tri2d(tri2,[xyz2t(:,1),xyz2t(:,2)],[0 0]);
   if ~in2
     ctr = mean(xyz2t);
     xyz1t = xyz1t-repmat(ctr,npts(k-1),1);
     xyz2t = xyz2t-repmat(ctr,npts(k),1);
   end
%
   in1 = in_tri2d(tri1,[xyz1t(:,1),xyz1t(:,2)],[0 0]);
   if ~in1
     error([' *** ERROR in mk_tri4p2:  Not able to find common', ...
            ' point in slices ' int2str(k-1) ' and ' int2str(k) '!']);
   end
%
% Transform to Polar Coordinates
%
   [th1,r1] = cart2pol(xyz1t(:,1),xyz1t(:,2));
   [th2,r2] = cart2pol(xyz2t(:,1),xyz2t(:,2));
%
% Get Directions of Angles
%
   dang = diff(th1);
   idx = find(abs(dang)<pi);           % Avoid the jump in angles
   dang = sum(dang(idx));
   dir1 = sign(dang);
%
   dang = diff(th2);
   idx = find(abs(dang)<pi);           % Avoid the jump in angles
   dang = sum(dang(idx));
   dir2 = sign(dang);
%
% Rotate Angles and Find Nearest Angle in Previous Slice
%
   th2 = th2-th1(1);
%
   idn = find(th2<0);
   th2(idn) = th2(idn)+2*pi;
%
   [~,idmn] = min(th2.*th2);
%
% Reorder Points
%
   if dir1~=dir2
     if dir2==1
       ido = [idmn-1:-1:1 npts(k):-1:idmn];
       idoc = [idmn ido];              % Including extra point to close slice
     else
       ido = [idmn:-1:1 npts(k):-1:idmn+1];
       idoc = [ido idmn];              % Including extra point to close slice
     end
   else
     if dir2==1
       ido = [idmn:npts(k) 1:idmn-1];
       idoc = [ido idmn];              % Including extra point to close slice
     else
       ido = [idmn+1:npts(k) 1:idmn];
       idoc = [idmn ido];              % Including extra point to close slice
     end
   end
%
   th2o = th2(ido);     % Ordered angles
%
% Add Point to Close Slice
%
   th2o = [th2o; (th2o(npts(k))+th2o(1))/2];% Add point to close slice
%
% Normalize "th1"
%
   th1 = [th1; (th1(npts(k-1))+th1(1))/2];  % Add point to close slice
   th1 = th1-th1(1);
   idn = find(th1<0);
   th1(idn) = th1(idn)+2*pi;
%
% Check Angles Are in Correct Order
%
   if dir1>0
     th1s = sort(th1,1,'ascend');
     th2os = sort(th2o,1,'ascend');
   else
     th1(1) = 2*pi;     % Equals zero - keeps position of points
     th1s = sort(th1,1,'descend');
     th2os = sort(th2o,1,'descend');
   end
%
% Slice Separations
%
   slx(k-1) = diff(cntr)*nv;
%
% Delaunay Triangulation
% Should the angles be scaled by arclength?
%
   xt = [zeros(npts(k-1)+1,1); slx(k-1)*ones(npts(k)+1,1)];
   yt = [th1s; th2os];
   tril = delaunay(xt,yt);
%
   if iplt
%    if ~iplt
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
% Check Normals
%
   xyz = [xyz1; xyz1(1,:); xyz2(idoc,:)];
   [nx,ny,nz,xc,yc,zc] = tri_norm(tril,xyz);
   ntril = size(tril,1);
   nv = [nx ny nz];
   xyzc = [xc yc zc]-repmat(ctr,ntril,1);
   irev = false(ntril,1);
   for l = 1:ntril
      irev(l) = xyzc(l,:)*nv(l,:)'<0;
   end
   if nnz(irev)>ntril/2 % If most are reversed, then all should be reversed
     tril = tril(:,[1 3 2]);
   end
%
% Improve Mesh
%
   tril = tri_fix2(tril,xyz);
%
% Global Node IDs
%
   nid1 = [n(k-1)+1:n(k) n(k-1)+1];
   nid2 = n(k)+1:n(k+1);
   nid2 = nid2(idoc);
   nid = [nid1 nid2];
   tri = [tri; nid(tril)];
%
% Create Triangles for Top and Bottom
%
   if k==2
     xyze = zeros(2,3);
     nidt = n(nslice+1)+1;             % Node ID for top center point
     xyze(1,:) = mean(dat{1});         % Coordinates for center point
%
     nptt = npts(k-1);
     trit = [repmat(nidt,nptt-1,1) [1:nptt-1]' [2:nptt]'];
     trit = [trit; [nidt nptt 1]];
%
     xyzt = zeros(n(nslice+1)+2,3);
     xyzt(1:nptt,:) = xyz1;
     xyzt(n(nslice+1)+1,:) = xyze(1,:);
     [nx,ny,nz,xc,yc,zc] = tri_norm(trit,xyzt);  % Check normals
     nv = [nx ny nz];
     xyzc = -diff(cntr);
     ntt = size(trit,1);
     irev = false(ntt,1);
     for l = 1:ntt
        irev(l) = xyzc*nv(l,:)'<0;
     end
     if nnz(irev)>ntt/2 % If most are reversed, then all should be reversed
       trit = trit(:,[1 3 2]);
     end
%
     tri = [tri; trit];
   end
%
   if k==nslice
     nidb = n(nslice+1)+2;             % Node ID for bottom center point
     xyze(2,:) = mean(dat{nslice});    % Coordinates for center point
%
     nptb = npts(k);
     trib = [repmat(nidb,nptb-1,1) [1:nptb-1]' [2:nptb]'];
     trib = [trib; [nidb nptb 1]];
%
     xyzb = zeros(n(nslice+1)+2,3);
     xyzb(1:nptb,:) = xyz2;
     xyzb(n(nslice+1)+2,:) = xyze(2,:);
     [nx,ny,nz,xc,yc,zc] = tri_norm(trib,xyzb);  % Check normals
     nv = [nx ny nz];
     xyzc = diff(cntr);
     ntb = size(trib,1);
     irev = false(ntb,1);
     for l = 1:ntb
        irev(l) = xyzc*nv(l,:)'<0;
     end
     if nnz(irev)>ntb/2 % If most are reversed, then all should be reversed
       trib = trib(:,[1 3 2]);
     end
%
     nid = n(k)+1:n(k+1);
     trib(:,2:3) = nid(trib(:,2:3));   % Global node IDs
     tri = [tri; trib];
   end
%
end
%
% Improve Mesh
%
xyz = cell2mat(dat);
xyz = [xyz; xyze];
tri = tri_fix2(tri,xyz);
%
nt = size(tri,1);
%
return