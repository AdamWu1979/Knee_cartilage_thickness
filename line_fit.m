function [pxyz,line_vec,v,score,pexp,res,sse] = line_fit(x,y,z,iplt);
%LINE_FIT  Fits a line to X, Y and Z point coordinate data using SVD.
%
%          [PXYZ,LINE_VEC] = LINE_FIT(X,Y,Z) given the X, Y and Z
%          coordinates of a set of points, calculates a least squares
%          fit of a line.  A point on the line, PXYZ, and the line
%          vector, LINE_VEC, give the parametric equation of the best
%          fit line by orthogonal regression.
%
%          [PXYZ,LINE_VEC,V,SCORE,PEXP,RES,SSE] = LINE_FIT(X,Y,Z)
%          returns the rotation matrix, V, PCA scores, SCORE, percent
%          of variance explained by the three orthogonal directions of
%          the plane, PEXP, the residuals (difference between the data
%          and fitted plane), RES, and the sum of squared errors, SSE.
%
%          NOTES:  1.  Must have at least three (3) points.
%
%                  2.  The SVD is used to do a principal component
%                  analysis (PCA) to do an orthogonal regression (total
%                  least squares) fit of the line.
%
%                  3.  Based on the demonstration of orthogonal
%                  regression using Matlab Statistics Toolbox.  See:
%                  http://www.mathworks.com/products/statistics/
%                  demos.html?file=/products/demos/shipping/stats/
%                  orthoregdemo.html
%
%          02-July-2010 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<4)
  iplt = 0;
end
%
if (nargin<3)
  error([' *** ERROR in LINE_FIT:  The X, Y and Z coordinates of ', ...
         'the points to be fit are required as inputs!']);
end
%
% Get Data Matrix and Check Number of Points
%
xyz = [x(:) y(:) z(:)];
npts = size(xyz,1);
if npts<3
  error(' *** ERROR in LINE_FIT:  Not enough data points!');
end
%
% Center Data
%
pxyz = mean(xyz);       % Point on the fitted line
xyz = xyz-repmat(pxyz,npts,1);         % Center data
%    
% Fit Line
%
[u,s,v] = svd(xyz);     % Number of datapoints x 3 plane parameters
%
line_vec = v(:,1);      % Line vector
%
% Additional Output
%
score = u*s;
%
if nargout>2
  rts = diag(s);
  pexp = 100*rts./sum(rts);
%
  res = xyz-score(:,1)*line_vec';
%
  sse = sum(sum(res.^2,2));
end
%
% Visual Check
%
if iplt
  figure;
  plot3(x,y,z,'bo','LineWidth',1);
  hold on;
  t = [min(score(:,1))-0.2; max(score(:,1))+0.2];
  line_fit = repmat(pxyz,2,1)+t*line_vec';
  plot3(line_fit(:,1),line_fit(:,2),line_fit(:,3),'k-','LineWidth',2);
  xyzf = repmat(pxyz,npts,1)+score(:,1)*line_vec';
  xyz = xyz+repmat(pxyz,npts,1);
  plot3([xyz(:,1) xyzf(:,1)]',[xyz(:,2) xyzf(:,2)]', ...
        [xyz(:,3) xyzf(:,3)]','r.-','LineWidth',1);
end
%
return