function [ip,il,ierr] = tsect4(v1,v2,v3,lp,lv,tol)
%TSECT4 Finds the intersection of a triangle and a line.
%
%       [IP,IL,IERR] = TSECT4(V1,V2,V3,LP,LV) finds the intersection of
%       a plane triangle defined by the vertices V1, V2 and V3 and a
%       line defined by a point LP and a direction vector LV.
%
%       TSECT4 returns the intersection point IP.  IL is set to true
%       (1) if an intersection was found.  IERR is set to true (1) if
%       the line lies (or is close to lying) in the plane of the
%       triangle.
%
%       [IP,IL,IERR] = TSECT4(V1,V2,V3,LP,LV,TOL) checks the
%       determinant to make sure it is not within tolerance TOL of zero
%       indicating the line is in (or close to) the plane of the
%       triangle.  Default tolerance is 1e-8.
%
%       NOTES:  1.  Based on the algorithm described in:
%               Tomas Moller and Ben Trumbore:  Fast, minimum storage
%               ray/triangle intersection.  Journal of Graphics Tools
%               2(1):21-28, 1997.
%
%               2.  Must have the M-file xprod.m in the current path or
%               directory.
%
%               3.  See also:
% http://www.mathworks.com/matlabcentral/fileexchange/25058-raytriangle-intersection
%
%               4.  Assumes the line is infinitely long.
%
%       10-Jul-2013 * Mack Gardner-Morse

%#######################################################################
%
% Check for Inputs
%
if (nargin<5)
  error(' *** Error in TSECT4:  Five input arguments are required!');
end
%
if (nargin<6)||isempty(tol)
  tol = 1e-8;
end
%
% Get Column Vectors
%
v1 = v1(:);
v2 = v2(:);
v3 = v3(:);
lp = lp(:);
lv = lv(:);
%
% Initial Values
%
il = false;
ierr = false;
ip = [];
%
% Calculate Edges
%
e1 = v2-v1;             % Edge 1
e2 = v3-v1;             % Edge 2
%
% Calculate Determinant
%
p = xprod(lv,e2);
d = p*e1;
%
% Check Determinant
%
if abs(d)<tol
  ierr = true;
  return
end
%
% Calculate Translation Vector
%
tv = lp-v1;             % Translation vector
%
% Calculate u and Test Bounds
%
u = p*tv/d;
if (u<0)||(u>1)
  return
end
%
% Calculate v and Test Bounds
%
q = xprod(tv,e1);
v = q*lv/d;
if (v<0)||(u+v>1)
  return
end
%
% Get Intersection Point
%
ip = (1-u-v)*v1+u*v2+v*v3;
il = true;
%
return