function xprd = xprod(v1,v2)
%XPROD  Calculates the cross product of two vectors.
%       XPRD = XPROD(V1,V2) Calculates the cross product of vectors V1
%       and V2.
%
%       NOTES:  1.  Returns the cross product as a row vector.
%
%       10-Jul-2013 * Mack Gardner-Morse

%#######################################################################
%
% Calculate the Cross Product
%
xprd = [v1(2).*v2(3)-v1(3).*v2(2), v1(3).*v2(1)-v1(1).*v2(3), ...
        v1(1).*v2(2)-v1(2).*v2(1) ];
%
return