function roi = rd_roi3(filenam,ipx);
%RD_ROI3  Reads an OSIRIX ROI CSV file.
%
%         ROI = RD_ROI3(FILENAM) reads the OSIRIX ROI CSV file, FILENAM,
%         and returns the structure, ROI, with the names of the ROI in
%         field "name" and the X, Y and Z data points for each slice is
%         in the columns of cell arrays in field "data". The X, Y and Z
%         data is in ordered triplets.  Each triplet represent a point
%         from a slice in the ROI.  Each ROI is returned in a row in
%         the structure.  
%
%         ROI = RD_ROI3(FILENAM,IPX) if IPX is true (nonzero) the X and
%         Y pixel data points for each slice is in the columns of cell
%         arrays in field "data" in the structure ROI.
%
%         NOTES:  1.  The data was collected from OSIRIX using the
%                 pencil tool.
%
%                 2.  If the ROI name doe not start with an alphabetic
%                 letter, the letter "A" is prepended to the ROI name.
%
%         16-June-2010 * Mack Gardner-Morse
%
%         03-May-2012 * Mack Gardner-Morse * Added flag and code to
%         return X and Y pixel data instead of X, Y and Z coordinate
%         data.
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<1)
  error(' *** ERROR in RD_ROI3:  An input file name is required!');
end
%
if (nargin<2)
  ipx = false;
end
%
if isempty(ipx)
  ipx = false;
end
%
% Open File and Read First Two Lines
%
fid = fopen(filenam,'rt');
lin = fgetl(fid);       % Read first line of headers
idx = strfind(lin,',');
idx = [idx length(lin)+1];
nwfm = strcmp(lin(idx(11)+1:idx(12)-1),'Length');
lin = fgetl(fid);       % Read first line of data
%
% Keep Reading Data Until End of File
%
while lin~=-1
     idx = strfind(lin,',');
     idx = [idx length(lin)+1];
     rnam = lin(idx(7)+1:idx(8)-1);    % Name of ROI
%
% Strip Any Quotes from ROI Name
%
     idq = rnam==''''|rnam=='"';       % Both single and double quotes
     idq = ~idq;
     rnam = rnam(idq);
%
% Check for Hyphens, Percents and Spaces
%
     htrap = findstr(rnam,'-');        % Trap for hyphens
     if ~isempty(htrap)
       rnam(htrap) = '_';              % Replace hyphens w/ underscore
     end
     htrap = findstr(rnam,'%');        % Trap for percents
     if ~isempty(htrap)
       rnam(htrap) = 'p';              % Replace percents w/ p's
     end
     htrap = findstr(rnam,' ');        % Trap for spaces
     if ~isempty(htrap)
       rnam(htrap) = '_';              % Replace spaces w/ underscore
     end
%
% Check that First Character is a Letter
%
     if ~isletter(rnam(1))
       rnam = ['a' rnam];              % Prepend an "a" if not a character
     end
%
% Number of Data Points for this ROI
%
     if nwfm
         npts = eval(lin(idx(14)+1:idx(15)-1));
         idp = 15+[0:5:(npts-1)*5]';
     else
         npts = eval(lin(idx(13)+1:idx(14)-1));
         idp = 14+[0:5:(npts-1)*5]';
     end
     if ipx
       mat = zeros(npts,2);
     else
       mat = zeros(npts,3);
     end
%
% Get Point Data
%
     for k = 1:npts;
        if ipx
          mat(k,:) = eval(['[' lin([idx(idp(k)+3)+1:idx(idp(k)+5)-1]) ']']);
        else
          mat(k,:) = eval(['[' lin([idx(idp(k))+1:idx(idp(k)+3)-1]) ']']);
        end
     end
%
% Save Data for Each ROI
%
     if exist(rnam,'var');
       eval(['nslice = size(' rnam ',2);']);
       eval([rnam '{' int2str(nslice+1) '} = mat;']);
     else
       eval([rnam '{1} = mat;']);
     end
     lin = fgetl(fid);
end
%
% Close File
%
fclose(fid);
%
% Clear Workspace
%
clear ans fid filenam htrap idq lin k idx ipx nslice rnam npts idp mat nwfm;
%
% Get ROI
%
var = whos;
rnam = {var.name}';     % Names of ROI
nvar = size(rnam,1);    % Number of ROI
%
% Get Data into a Cell
%
data = cell(nvar,1);
for k = 1:nvar
%   eval(['roi(' int2str(k) ').data = eval( char( rnam(' int2str(k) ')));']);
  eval(['data{' int2str(k) '} = eval( char( rnam(' int2str(k) ')));']);
end
%
% Put Information into a Structure
%
roi = struct('name',rnam,'data',data);
%
return