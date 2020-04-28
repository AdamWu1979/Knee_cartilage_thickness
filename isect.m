function [ip,ierr,t] = isect(pp,pv,lp,lv,tol)
%ISECT  Finds the intersection of a plane and a line.
%       [IP,IERR] = ISECT(PP,PV,LP,LV) finds the intersection of a
%       plane defined by a point (PP) and a normal vector (PV) and a
%       line defined by a point (LP) and a direction vector (LV).
%       ISECT returns the intersection point (IP).  IERR is set to
%       one (1) if there was an error.
%
%       IP = ISECT(PP,PV,LP,LV,TOL) checks the intersection is within
%       a tolerance of the plane.  Default tolerance is 1e-8.
%
%       [IP,IERR,T] = ISECT(PP,PV,LP,LV,TOL) returns the distance (T)
%       along the vector (LV) to the intersection with the plane.
%
%       29-Mar-96 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  error('ISECT requires four input arguments.');
end
%
if (nargin<5)
  tol = 1e-8;
end
%
% Error Flag
%
ierr = 0;
%
% Check Vectors
%
pp = pp(:);
pv = pv(:);
lp = lp(:);
lv = lv(:);
%
[n1 l1] = size(pp);
[n2 l2] = size(pv);
[n3 l3] = size(lp);
[n4 l4] = size(lv);
%
if ((l1~=1)|(l2~=1)|(l3~=1)|(l4~=1))
  error('ISECT only works with vectors.')
end
%
% Check that the Inputs have Three Rows
%
if ((n1~=3)|(n2~=3)|(n3~=3)|(n4~=3))
  error('Point and vector inputs must be of length three (3).')
end
%
% Find Intersection
%
d = pv'*lv;
if (abs(d)>10*eps)
  t = lp-pp;
  t = -pv'*t;
  t = t/d;
%
  ip = t*lv+lp;
%
% Check Intersection
%
  zo = pv'*(ip-pp);
  if (abs(zo)>tol)
    ierr = 1;
    warning(' *** WARNING in ISECT:  Intersection not in plane.')
  end
else
  ierr = 1;
  ip = [];
  t = [];
  warning(' *** WARNING in ISECT:  Divide by small number.');
end
%
return