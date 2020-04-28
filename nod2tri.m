function [ir,nt] = nod2tri(nodlst,tric,ncn);
%NOD2TRI  Finds the triangles connected to a vector of nodes.
%
%         IR = NOD2TRI(NODLST,TRIC) returns the row index to triangles
%         connected to the nodes in the vector, NODLST, that are in the
%         three (3) column triangle connectivity matrix, TRIC.
%
%         [IR,NT] = NOD2TRI(NODLST,TRIC) returns the number of
%         triangles, NT, connected to the nodes in NODLST.
%
%         IR = NOD2TRI(NODLST,TRIC,NCN) returns the row index to
%         triangles connected to the nodes in NODLST with greater
%         than NCN nodes from NODLST in the triangle.  NCN must be
%         between zero (0) and two (2) since there are only three
%         nodes in a triangle.  The default is zero (0) (at least
%         one node from the triangle is in the input list).
%
%         NOTES:  1.  Note that the number of connected nodes from
%                 NODLST must be greater than the number of connected
%                 nodes, NCN, in a triangle.
%
%         27-Sep-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check Inputs
%
if (nargin<2)
  error(' *** ERROR in NOD2TRI:  Not enough input arguments!');
end
%
if (nargin<3)
  ncn = 0;
end
%
if size(tric,2)~=3
  error([' *** ERROR in NOD2TRI:  Triangle connectivity matrix,', ...
         ' TRIC, must have three (3) columns!']);
end
%
if ncn<0|ncn>2
  error([' *** ERROR in NOD2TRI:  Number of connected nodes, NCN,' ...
         ' must be between zero (0) and two (2)!']);
end
%
% Get Sizes of Input Data
%
nodlst = nodlst(:);
nn = size(nodlst,1);
%
nt = size(tric,1);
%
% Find Connected Triangles
%
nnc = false(nt,3);
for k = 1:nn
   nnc = nnc|tric==nodlst(k);
end
ir = find(sum(nnc,2)>ncn);
nt = size(ir,1);
%
return