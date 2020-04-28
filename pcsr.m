function [cgs,rotmat,xyzlt2,veclt2,vols,tris,xyzs,inertia_cgs] = ...
          pcsr(leg,dats,datx,iplt);
%PCSR     Determines the patella coordinate system (PCS) using the
%         sagittal and axial digitized patella bone MRI using the
%         method of Rainbow et al.
%
%         [CGS,ROTMAT] = PCSR(LEG,DATS,DATX) given a logical scalar,
%         LEG (1 [true] for right and 0 [false] for left patellas), a
%         cell array containing three (3) columns matrices with slice
%         coordinate point data for the sagittal digitization, DATS,
%         and a cell array containing three (3) columns matrices with
%         slice coordinate point data for the axial digitization, DATX,
%         returns a three (3) column row vector with the origin (center
%         of gravity) and a 3x3 rotation matrix, ROTMAT, from the MRI
%         to patella coordinate system (PCS).
%
%         [CGS,ROTMAT,XYZLT2,VECLT2] = PCSR(LEG,DATS,DATX) returns the
%         line fitted to the ridge points in the patella coordinate
%         system as a center point, XYZLT2, and a vector, VECLT2.
%
%         [CGS,ROTMAT,XYZLT2,VECLT2,VOLS] = PCSR(LEG,DATS,DATX) returns
%         the volume of the sagittal digitization in VOLS.
%
%         [CGS,ROTMAT,XYZLT2,VECLT2,VOLS,TRIS,XYZS] = PCSR(LEG,DATS,
%         DATX) returns the triangle connectivity matrix, TRIS, and the
%         X, Y and Z coordinate matrix, XYZS, in the patella coordinate
%         system.
%
%         [CGS,ROTMAT,XYZLT2,VECLT2,VOLS,TRIS,XYZS,INERTIA_CGS] = PCSR(
%         LEG,DATS,DATX) returns the inertia matrix, INERTIA_CGS, which
%         is the inertia about the center of gravity in the MRI
%         coordinate system.
%
%         NOTES:  1.  Each slice coordinate data matrix must correspond
%                 to one index into the MRI data cell arrays (DATS and
%                 DATX).  Assumes the slice data is in sequential order.
%
%                 2.  The M-files coord_tf.m, in_tri2d.m, inert_tri.m,
%                 isect.m, line_fit.m, mk_tri4p2.m, plane_fit.m,
%                 rd_roi3.m, sl_info.m, tri_fix2.m and tri_norm.m must
%                 be in the current path or directory.
%
%                 3.  The methodology is based on the method outlined
%                 in:  Rainbow MJ, Miranda DL, Cheung RT, Schwartz JB,
%                 Crisco JJ, Davis IS, Fleming BC.  Automatic
%                 determination of an anatomical coordinate system for
%                 a three-dimensional model of the human patella.
%                 J Biomech 46(12):2093-6, 2013.  PMID:  23791087
%                 doi: 10.1016/j.jbiomech.2013.05.024.
%
%                 4.  The Simulink 3D Animation toolbox is required.
%
%         14-Jul-2017 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR in PCSR:  Not enough input data!');
end
%
if (nargin<4)
  iplt = false;
end
%
if ~isscalar(leg)
  error(' *** ERROR in PCSR:  LEG variable is not scalar!');
end
%
leg = logical(leg(1));
%
if ~iscell(dats)||~iscell(datx)
  error(' *** ERROR in PCSR:  MRI data is not in cell arrays!');
end
%
% Algorithm Parameters
% (Make Input Variables?)
%
npts = 3;               % Number of points to define ridge in slice
% rchk2 = 36;             % Find mean slice points more than 6 mm from the fitted line
rchk2 = 25;             % Find mean slice points more than 5 mm from the fitted line
%
% Generate Mesh for Sagittal Data
%
[tris,xyzs] = mk_tri4p2(dats);
%
% Get Inertias and Eigenvectors
%
[vols,cgs,~,inertia_cgs] = inert_tri(xyzs,tris);
%
[evecs,~] = eig(inertia_cgs);
%
% Get AP Vector
%
[~,aplocs] = max(abs(evecs(2,:)));     % Get vector in MRI Y-direction
apvecs = evecs(:,aplocs)';
%
if apvecs(2)>0
  apvecs = -apvecs;
end
%
% Get Rotation Vector and Angle and First Rotation Matrix
%
if leg
  rotang1 = vrrotvec(apvecs,[-1 0 0]);
else
  rotang1 = vrrotvec(apvecs,[1 0 0]);
end
%
rotmat1 = vrrotvec2mat(rotang1)';
%
% Plot Sagittal Data and Eigenvector
%
if iplt
  figure;
  orient tall;
  plt_datsl(dats,'b.-',0.5);
  hold on;
  sc = 20;              % Length of vector
  hx1 = quiver3(cgs(1),cgs(2),cgs(3),sc*apvecs(1),sc*apvecs(2), ...
                sc*apvecs(3),0,'r-');  % Centroid and eigenvector
  set(hx1,'LineWidth',3);
  plot3(cgs(1),cgs(2),cgs(3),'rs','LineWidth',2, ...
        'MarkerFaceColor','r','MarkerSize',10);
  xlabel('MRI X (mm)','FontSize',12,'FontWeight','bold');
  ylabel('MRI Y (mm)','FontSize',12,'FontWeight','bold');
  zlabel('MRI Z (mm)','FontSize',12,'FontWeight','bold');
  view(3);
  axis equal;
  title('MRI Coordinates','FontSize',16,'FontWeight','bold');
%   title({['Knee ' kid];['MRI Coordinates']},'Interpreter','none', ...
%         'FontSize',16,'FontWeight','bold');
%
%   if isave
%     psnam = ['pcsr_' kid '.ps'];
%     print('-dpsc2',psnam);
%   end
%
  datst1 = coord_tf(cgs,rotmat1,dats); % Transform data
  apvecst1 = apvecs*rotmat1;
%
  figure;
  orient tall;
  plt_datsl(datst1,'b.-',0.5);
  hold on;
  hx2 = quiver3(0,0,0,sc*apvecst1(1),sc*apvecst1(2), ...
                sc*apvecst1(3),0,'rs-','filled');  % Centroid and eigenvector
  set(hx2,'LineWidth',3);
  plot3(0,0,0,'rs','LineWidth',2,'MarkerFaceColor','r', ...
        'MarkerSize',10);
  xlabel('X (mm)','FontSize',12,'FontWeight','bold');
  ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
  zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
  view(40,25);
  axis equal;
  title('First Rotation','FontSize',16,'FontWeight','bold');
%   title({['Knee ' kid];['First Rotation']},'Interpreter','none', ...
%         'FontSize',16,'FontWeight','bold');
%
%   if isave
%     print('-dpsc2','-append',psnam);
%   end
%
end
%
% Transform Axial Data to Find the Patella Ridge
%
datxt1 = coord_tf(cgs,rotmat1,datx);
%
% Find Patella Ridge
%
nsl = size(datxt1,1);   % Number of slices
%
xyzr = NaN(npts*nsl,3); % Initial ridge point array
%
for l = 1:nsl
   idl = npts*l;
   idl = idl-npts+1:idl;
%
   xyzsl = datxt1{l};
%
   if size(xyzsl,1)>npts
     if leg
       [~,ids] = sort(xyzsl(:,1),'descend');
     else
       [~,ids] = sort(xyzsl(:,1),'ascend');
     end
     xyzr(idl,:) = xyzsl(ids(1:npts),:);         % Ridge points
   end
end
%
idg = ~isnan(xyzr(:,1));               % Check for missing values
xyzr = xyzr(idg,:);
%
% Find Change in Ridge Heights (Not Used)
%
if leg
  [~,msl] = max(xyzr(:,1));
else
  [~,msl] = min(xyzr(:,1));
end
msls = npts*ceil(msl/npts);
ide = 1:msls;        % Index to good ridge points
%
% Plot Ridge Points
%
if iplt
  figure;
  orient tall;
  plt_datsl(datxt1,'b.-',0.5);
  hold on;
  plot3(xyzr(msls+1:end,1),xyzr(msls+1:end,2),xyzr(msls+1:end,3), ...
        'rs','MarkerSize',8,'LineWidth',1);
  plot3(xyzr(1:msls,1),xyzr(1:msls,2),xyzr(1:msls,3),'ks', ...
        'MarkerSize',8,'LineWidth',1);
  xlabel('X (mm)','FontSize',12,'FontWeight','bold');
  ylabel('Y (mm)','FontSize',12,'FontWeight','bold');
  zlabel('Z (mm)','FontSize',12,'FontWeight','bold');
  view(3);
  axis equal;
end
%
% Fit Line to Ridge Points
%
% xyzr = xyzr(ide,:);
[xyzl,vecl,~,~,~,res] = line_fit(xyzr(:,1),xyzr(:,2),xyzr(:,3));
%
% Check Fit
%
d2res = sum(res.*res,2);               % Squared distance to line
npte = size(xyzr,1);    % Number of good points
nslg = npte/npts;       % Number of good slices
d2res = reshape(d2res,npts,nslg);
d2resmn = mean(d2res)'; % Mean squared distance to fitted line for each slice
idg = true(npts,nslg);
% idb = find(d2resmn>49); % Find mean slice points more than 7 mm from the fitted line
idb = find(d2resmn>rchk2);             % Find mean slice points more than sqrt (RCHK2) mm from the fitted line
while ~isempty(idb)
     idg(:,idb) = false;
     idg = idg(:);
     xyzr = xyzr(idg,:);
     [xyzl,vecl,~,~,~,res] = line_fit(xyzr(:,1),xyzr(:,2),xyzr(:,3));     % New fit
     d2res = sum(res.*res,2);          % Squared distance to line
     npte = size(xyzr,1);              % Number of good points
     nslg = npte/npts;                 % Number of good slices
     d2res = reshape(d2res,npts,nslg);
     d2resmn = mean(d2res)';           % Mean squared distance to fitted line for each slice
     idg = true(npts,nslg);
     idb = find(d2resmn>rchk2);        % Find mean slice points more than sqrt (RCHK2) mm from the fitted line
end
%
% Plot Final Ridge Points and Fitted Line
%
if iplt
  plot3(xyzr(:,1),xyzr(:,2),xyzr(:,3),'go','MarkerSize',8, ...
        'LineWidth',1);
  rdist = mean(xyzr(1:npts,:))-mean(xyzr(npte-npts+1:npte,:));
  rdist = 1.2*norm(rdist)/2;           % Half distance between ends of ridge points
  plot3([xyzl(1)-rdist*vecl(1); xyzl(1)+rdist*vecl(1)], ...
        [xyzl(2)-rdist*vecl(2); xyzl(2)+rdist*vecl(2)], ...
        [xyzl(3)-rdist*vecl(3); xyzl(3)+rdist*vecl(3)], ...
        'r-','LineWidth',2);
  view(-40,10);
  axis equal;
  title('Ridge Points','FontSize',16,'FontWeight','bold');
%   title({['Knee ' kid];['Ridge Points']},'Interpreter','none', ...
%         'FontSize',16,'FontWeight','bold');
%   if isave
%     print('-dpsc2','-append',psnam);
%   end
end
%
% Get Rotation Vector and Angle and Second Rotation Matrix
%
if vecl(3)<0
  vecl = -vecl;
end
%
xax = [1 0 0];
yax = cross(vecl,xax);
yax = yax/norm(yax);
zax = cross(xax,yax);
zax = zax/norm(zax);
xax = cross(yax,zax);
xax = xax/norm(xax);
rotmat2 = [xax; yax; zax]';
%
% Final Rotation Matrix
%
rotmat = rotmat1*rotmat2;
%
if nargout>2
  veclt2 = (vecl'*rotmat2)';
  xyzlt2 = xyzl*rotmat2;
  nps = size(xyzs,1);
  xyzs = (xyzs-repmat(cgs,nps,1))*rotmat;        % Transform coordinates
end
%
if iplt
%
  datst2 = coord_tf(cgs,rotmat,dats);  % Transform data
  apvecst2 = apvecs*rotmat;
%
  datxt2 = coord_tf(zeros(1,3),rotmat2,datxt1);
  veclt2 = (vecl'*rotmat2)';
  xyzlt2 = xyzl*rotmat2;
%
  figure;
  orient tall;
  plt_datsl(datxt2,'b.-',0.5);
  hold on;
  plt_datsl(datst2,'k.-',0.5);
%
  hx3 = quiver3(0,0,0,sc*apvecst2(1),sc*apvecst2(2), ...
                sc*apvecst2(3),0,'rs-','filled');     % Centroid and eigenvector
  set(hx3,'LineWidth',3);
  plot3(0,0,0,'rs','LineWidth',2,'MarkerFaceColor','r', ...
        'MarkerSize',10);
%
  plot3([xyzlt2(1)-rdist*veclt2(1); xyzlt2(1)+rdist*veclt2(1)], ...
        [xyzlt2(2)-rdist*veclt2(2); xyzlt2(2)+rdist*veclt2(2)], ...
        [xyzlt2(3)-rdist*veclt2(3); xyzlt2(3)+rdist*veclt2(3)], ...
        'r-','LineWidth',3);
%
  xlabel('Patella X (mm)','FontSize',12,'FontWeight','bold');
  ylabel('Patella Y (mm)','FontSize',12,'FontWeight','bold');
  zlabel('Patella Z (mm)','FontSize',12,'FontWeight','bold');
  view(3);
  axis equal;
  title('Patella Coordinates','FontSize',16,'FontWeight','bold');
%   title({['Knee ' kid];['Patella Coordinates']},'Interpreter', ...
%         'none','FontSize',16,'FontWeight','bold');
%
%   if isave
%     print('-dpsc2','-append',psnam);
%   end
end
%
return