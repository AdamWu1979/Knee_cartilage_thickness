function [vol,cg,inertia,inertia_cg] = inert_tri(xyz,tri);
%INERT_TRI Calculates the three-dimensional inertia of an enclosed
%          triangular meshed surface assuming an uniform density for
%          the enclosed volume.
%
%          [VOL,CG,INERTIA,INERTIA_CG] = INERT_TRI(XYZ,TRI) - given the
%          three column matrix of X, Y and Z coordinate, XYZ, and a
%          three column triangle connectivity matrix, TRI, returns the
%          enclosed volume, VOL, center of gravity (CG), moment of
%          inertia matrix about the coordinate system origin, INERTIA,
%          and momemt of inertia matrix about the center of gravity,
%          INERTIA_CG.
%
%          NOTES:  1.  Surface integral algorithm based on David
%                  Eberly's Geometric Tools, Redmond WA 98052.  See
%                  links for method details:
%                  https://www.geometrictools.com/
% http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
% http://www.geometrictools.com/GTEngine/Include/GtePolyhedralMassProperties.h
%
%          17-Sep-2015 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<2)
  error(' *** ERROR in INERT_TRI:  Two input variables are required!');
end
%
% Check Inputs
%
ncc = size(xyz,2);
[nt,nct] = size(tri);
%
if ncc~=3
  error([' *** ERROR in INERT_TRI:  Input coordinate matrix must', ...
         ' have three (3) columns!']);
end
%
if nct~=3
  error([' *** ERROR in INERT_TRI:  Input triangle connectivity', ...
         ' matrix must have three (3) columns!']);
end
%
% 3D Enclosed Triangular Meshed Surface Method
% See:  http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
% http://www.geometrictools.com/GTEngine/Include/GtePolyhedralMassProperties.h
%
o = ones(3,1);
rc = [20; 5*o; 2*o; o]/120;
ri = zeros(10,1);
%
for k = 1:nt
%
   v1 = xyz(tri(k,1),:);
   v2 = xyz(tri(k,2),:);
   v3 = xyz(tri(k,3),:);
   s1 = v2-v1;
   s2 = v3-v1;
   nv = cross(s1,s2);
%
   tmp1 = v1+v2;
   f1 = tmp1+v3;
   tmp2 = v1.*v1;
   tmp3 = tmp2+v2.*tmp1;
   f2 = tmp3+v3.*f1;
   f3 = v1.*tmp2+v2.*tmp3+v3.*f2;
   g1 = f2+v1.*(f1+v1);
   g2 = f2+v2.*(f1+v2);
   g3 = f2+v3.*(f1+v3);
%
   ri(1) = ri(1)+nv(1)*f1(1);
   ri(2:4) = ri(2:4)+(nv.*f2)';
   ri(5:7) = ri(5:7)+(nv.*f3)';
   ri(8:9) = ri(8:9)+(nv(1:2).*(v1(2:3).*g1(1:2)+v2(2:3).*g2(1:2)+ ...
             v3(2:3).*g3(1:2)))';
   ri(10) = ri(10)+(nv(3).*(v1(1).*g1(3)+v2(1).*g2(3)+v3(1).*g3(3)));
%  
end
%
ri = rc.*ri;
%
vol = ri(1);            % Volume
cg = ri(2:4)'./vol;     % Center of gravity (CG)
inertia = [ri(6)+ri(7) -ri(8) -ri(10); -ri(8) ri(5)+ri(7) -ri(9); ...
          -ri(10) -ri(9) ri(5)+ri(6)]; % Inertia
%
iicg = cg*cg'*eye(3)-cg'*cg;           % Parallel axis theorem
%
inertia_cg = inertia-vol*iicg;         % Inertia about CG
%
return