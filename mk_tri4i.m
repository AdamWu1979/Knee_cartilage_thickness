function [tri,nt,slx,ierr] = mk_tri4i(dat,imprv,lcirc,iends,tol,iplt);
%MK_TRI4I Makes a triangular mesh by using the ordered slice data from
%        the digitized MRI using arclengths along the slices.  The
%        function includes options to close slices and enclose the
%        first and last slices.  The function also tries to improve the
%        triangle aspect ratios.
%
%        [TRI,NT] = MK_TRI4I(DAT) given a cell array containing three (3)
%        columns matrices with slice coordinate point data, DAT, returns
%        the three (3) column triangle connectivity matrix, TRI.  The
%        number of returned triangles, NT, may also be returned.
%
%        TRI = MK_TRI4I(DAT,IMPRV) if logical IMPRV is true (nonzero),
%        function tries to improve the triangle aspect ratios.  Use
%        cautiously as the algorithm may change the mesh topology.  By
%        default, the function tries to improve the triangle aspect
%        ratios.
%
%        [TRI,NT,SLX,IERR] = MK_TRI4I(DAT,IMPRV) returns the slice
%        separation distances in column vector SLX and IERR is true if
%        there was an error in trying to improve the triangle aspect
%        ratios of the final full mesh.
%
%        TRI = MK_TRI4I(DAT,IMPRV,LCIRC) if logical LCIRC is true,
%        will connect the first and last points in the slice data to
%        make a mesh with a closed bounadry.  By default, the function
%        does not connect the first and last points.
%
%        TRI = MK_TRI4I(DAT,IMPRV,LCIRC,IENDS) when logical LCIRC is
%        true, and integer IENDS is between 1 and 3 will mesh the ends
%        (first and/or last slice data).  The number of end meashes
%        depends on the value of IENDS (0 = none, 1 = top, 2 = bottom
%        and 3 = both).  The ends are not enclosed by default.
%
%        NOTES:  1.  Each slice coordinate data matrix must correspond
%                to one index into the cell array DAT.
%
%                2.  The coordinates should be ordered in the same
%                direction in every slice.  The dot product of the
%                directions of adjacent slices are used to check the
%                ordering direction and the ordering direction is
%                reversed if the dot product is negative.
%
%                3.  The arclength along each slice is used to determine
%                the triangulation.
%
%                4.  The algorithm that tries to improve triangle
%                aspect ratios may change the mesh topology.
%
%                5.  The triangle normals may not point outwards for
%                odd shaped slice boundaries.
%
%                6.  The end meshes assume that the majority of 
%                triangles are within the slice boundary.
%
%                7.  The M-file nod2tri.m, plane_fit.m, tri_fix2.m and
%                tri_norm.m must be in the current path or directory.
%
%        22-Jan-2016 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<6)||isempty(iplt)
  iplt = false;
end
%
if iplt
  hf = figure;
  orient tall;
end
%
if (nargin<5)||isempty(tol)
  tol = 0.1;            % Tolerance on dot product (projection of directions)
end
%
if (nargin<4)||isempty(iends)
  iends = 0;            % Number of end caps - 0 = none, 1 = top, 2 = bottom and 3 = both
end
%
if (nargin<3)||isempty(lcirc)
  lcirc = false;        % Slices are loops?
end
%
if (nargin<2)||isempty(imprv)
  imprv = true;         % Improve triangle aspect ratios
end
%
if (nargin<1)
  error(' *** ERROR in MK_TRI4I:  No input data!');
end
%
if lcirc
  iends = round(iends);
  if iends>3||iends<0
    iends = 0;
    warning(' *** WARNING in MK_TRI4I:  IENDS value not recognized!');
  end
end
%
% Get Arc Lengths
%
dat = dat(:);
nslice = size(dat,1);
slen = cell(nslice,1);
npts = zeros(nslice,1);
rpt1 = zeros(nslice,3);
rpt2 = zeros(nslice,3);
vec1 = zeros(nslice,3);
cntr = zeros(nslice,3);
nvec = zeros(3,nslice);
irev = false(nslice,1);
%
for k = 1:nslice
   xyz = dat{k};
   [cntr(k,:),nvec(:,k)] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));
   vec = xyz(end,:)-xyz(1,:);
   vec1(k,:) = vec./norm(vec);
%
% Check for Slices with a Reverse Digitization
%
   if k>1
     dotp = vec1(k-1,:)*vec1(k,:)';
     if dotp<tol
       irev(k) = true;
       xyz = flipud(xyz);
       vec = xyz(end,:)-xyz(1,:);
       vec1(k,:) = vec./norm(vec);
       dotp2 = vec1(k-1,:)*vec1(k,:)';
       if dotp2<dotp    % Revert back to original ordering
         warning([' *** WARNING in mk_tri4i:  Ordering of points', ...
                  ' in the slices may not be in the same direction!']);
         irev(k) = false;
         xyz = flipud(xyz);
         vec = xyz(end,:)-xyz(1,:);
         vec1(k,:) = vec./norm(vec);
       end
     else
       irev(k) = false;
     end
   end
   rpt1(k,:) = xyz(1,:);
   rpt2(k,:) = xyz(2,:);
   npts(k) = size(xyz,1);
   if lcirc
     xyz = [xyz; rpt1(k,:)];
   end
   dd = diff(xyz);
   dlen = sqrt(sum(dd.*dd,2));
   if irev(k)
     slen{k} = flipud([0; cumsum(dlen)]);
   else
     slen{k} = [0; cumsum(dlen)];
   end
end
%
n = [0; cumsum(npts)];
tri = [];
slx = zeros(nslice-1,1);
for k = 2:nslice
%
% Slice Separations and Offsets
%
   ds = rpt1(k,:)-rpt1(k-1,:);
   slx(k-1) = (cntr(k,:)-cntr(k-1,:))*nvec(:,k-1);
   if lcirc
     offst = ds*(rpt2(k-1,:)-rpt1(k-1,:))';
   else
     offst = ds*vec1(k-1,:)';
   end
%
% Delaunay Triangulation
%
   if lcirc
     xt = [zeros(npts(k-1)+1,1); slx(k-1)*ones(npts(k)+1,1)];
   else
     xt = [zeros(npts(k-1),1); slx(k-1)*ones(npts(k),1)];
   end
   yt = [slen{k-1}-offst; slen{k}];
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
% Get Slice Coordinates
%
   xyz1 = dat{k-1};
   xyz2 = dat{k};
   if lcirc
     xyz1 = [xyz1; xyz1(1,:)];
     xyz2 = [xyz2; xyz2(1,:)];
   end
   if irev(k-1)
     xyz1 = flipud(xyz1);
   end
   if irev(k)
     xyz2 = flipud(xyz2);
   end
%
   xyz = [xyz1; xyz2];
%
% Check Normals
%
   [nx,ny,nz,xc,yc,zc] = tri_norm(tril,xyz);
   ntril = size(tril,1);
   nv = [nx ny nz];
%
   ctr = mean(cntr(k-1:k,:));
   xyzc = [xc yc zc]-repmat(ctr,ntril,1);
   irevn = false(ntril,1);
   for l = 1:ntril
      irevn(l) = xyzc(l,:)*nv(l,:)'<0;
   end
   if nnz(irevn)>ntril/2               % If most are reversed, then all should be reversed
     tril = tril(:,[1 3 2]);
   end
%
% Improve Mesh
%
   if imprv
     tril = tri_fix2(tril,xyz);
   end
%
% Global Node IDs
%
   if lcirc
     nid = [n(k-1)+1:n(k) n(k-1)+1 n(k)+1:n(k+1) n(k)+1];
   else
     nid = n(k-1)+1:n(k+1);
   end
%
   tri = [tri; nid(tril)];
%
% Create Triangles for Top Slices
%
   if k==2&&lcirc&&(iends==1||iends==3)
%
     xyz = dat{k-1};
     xyz = xyz-repmat(cntr(k-1,:),npts(k-1),1);
     [~,~,~,r] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));
     xyzt = xyz*r;      % Rotate slice into XY plane
     ntid = [1:npts(k-1) 1]';
     trit = delaunay(xyzt(:,1),xyzt(:,2));       % Delaunay triangulation
     [~,~,~,xc,yc,zc] = tri_norm(trit,xyzt);     % Get triangle centers
     nt = size(trit,1);
     t = zeros(nt,1);   % Which side of boundary?
%
% Use Cross Product to Determine Side of Boundary
%
     for l = 2:npts(k-1)
        pt1 = xyzt(l-1,:);
        v1 = xyzt(l,:)-pt1;
        it = nod2tri(ntid(l-1:l),trit,1);
        nt = size(it,1);
        vt = [xc(it) yc(it) zc(it)]-repmat(pt1,nt,1);
        v2 = cross(repmat(v1,nt,1),vt);
        t(it) = sign(v2(:,3));         % Check Z coordinate
     end
%
     id0 = find(t==0);  % Check triangles with no sides on the boundary (only vertices on the boundary)
     n0 = size(id0,1);
     tri0 = [trit(id0,:) trit(id0,1)];
     dt = abs(diff(tri0,1,2));         % Look for closest two points along boundary
     [~,mi] = min(dt,[],2);
%
     for l = 1:n0
        it = id0(l);
        m = sort(tri0(l,mi(l):mi(l)+1));
        pt1 = xyzt(m(1),:);
        v1 = xyzt(m(2),:)-pt1;
        vt = [xc(it) yc(it) zc(it)]-pt1;
        v2 = cross(v1,vt);
        t(it) = sign(v2(:,3));         % Check Z coordinate
     end
%
     if sum(t)>0        % Assume most triangles are within boundary
       it = find(t>0);
     else
       it = find(t<0);
     end
     trit = trit(it,:);
%
     [nx,ny,nz] = tri_norm(trit,xyz);  % Check normals
     nv = [nx ny nz];
     xyzc = -diff(cntr(k-1:k,:));
     ntt = size(trit,1);
     irevn = nv*xyzc'<0;
     if nnz(irevn)>ntt/2               % If most are reversed, then all should be reversed
       trit = trit(:,[1 3 2]);
     end
%
     tri = [tri; trit];
   end
%
% Create Triangles for Bottom Slices
%
   if k==nslice&&lcirc&&(iends==2||iends==3)
%
     xyz = dat{k};
     xyz = xyz-repmat(cntr(k,:),npts(k),1);
     [~,nvpl,~,r] = plane_fit(xyz(:,1),xyz(:,2),xyz(:,3));
     xyzt = xyz*r;      % Rotate slice into XY plane
     ntid = [1:npts(k) 1]';
     trib = delaunay(xyzt(:,1),xyzt(:,2));       % Delaunay triangulation
     [~,~,~,xc,yc,zc] = tri_norm(trib,xyzt);     % Get triangle centers
     nt = size(trib,1);
     t = zeros(nt,1);   % Which side of boundary?
%
% Use Cross Product to Determine Side of Boundary
%
     for l = 2:npts(k)
        pt1 = xyzt(l-1,:);
        v1 = xyzt(l,:)-pt1;
        it = nod2tri(ntid(l-1:l),trib,1);
        nt = size(it,1);
        vt = [xc(it) yc(it) zc(it)]-repmat(pt1,nt,1);
        v2 = cross(repmat(v1,nt,1),vt);
        t(it) = sign(v2(:,3));         % Check Z coordinate
     end
%
     id0 = find(t==0);  % Check triangles with no sides on the boundary (only vertices on the boundary)
     n0 = size(id0,1);
     tri0 = [trib(id0,:) trib(id0,1)];
     dt = abs(diff(tri0,1,2));         % Look for closest two points along boundary
     [~,mi] = min(dt,[],2);
%
     for l = 1:n0
        it = id0(l);
        m = sort(tri0(l,mi(l):mi(l)+1));
        pt1 = xyzt(m(1),:);
        v1 = xyzt(m(2),:)-pt1;
        vt = [xc(it) yc(it) zc(it)]-pt1;
        v2 = cross(v1,vt);
        t(it) = sign(v2(:,3));         % Check Z coordinate
     end
%
     if sum(t)>0        % Assume most triangles are within boundary
       it = find(t>0);
     else
       it = find(t<0);
     end
     trib = trib(it,:);
%
     [nx,ny,nz] = tri_norm(trib,xyz);  % Check normals
     nv = [nx ny nz];
     xyzc = diff(cntr(k-1:k,:));
     ntb = size(trib,1);
     irevn = nv*xyzc'<0;
     if nnz(irevn)>ntb/2               % If most are reversed, then all should be reversed
       trib = trib(:,[1 3 2]);
     end
%
     nid = n(k)+1:n(k+1);
     trib = nid(trib);  % Global node IDs
     tri = [tri; trib];
   end
%
end
%
% Improve Mesh
%
if imprv
  xyz = cell2mat(dat);
  try
     tri = tri_fix2(tri,xyz);
     ierr = false;
  catch
     warning(' *** WARNING in MK_TRI4I:  Unable to improve mesh!');
     ierr = true;
  end
end
%
nt = size(tri,1);
%
return