function mat2vec(x,y,u,v,fname,varargin)
% MAT2VEC converts Matlab matrix of the specific
% oceonagrophic data set into VEC files, suitable 
% for the GRADPIV2 and TECPLOT software.
%
% Inputs: 
%       X,Y - two 2D matrices of X and Y coordinates
%       U,V - two 3D matrices of U,V velocities, 3rd dimension 
%               is the number of flow maps.
%       BASENAME - part of the VEC filename, that will construct
%               the series of files, BASENAME0001.VEC, BASENAME0002.VEC,...
%
% Outputs:
%        
%
% Example:
%           load ocean500.mat
%           mat2vec(x,z,u2,w2,'OceanData');
%
%
% See also 
%
% Author: Alex Liberzon
% Date: 19-Apr-2004
%
% Last modified:
% Version: 

% we need x,y as rows
x = x(:)';
y = y(:)';
u = u';
v = v';

r = length(unique(x));
c = length(unique(y));

for i = 1:size(u,1)
    fid = fopen(sprintf('%s%04d.vec',fname,i),'w');
    header = sprintf('TITLE="%s" VARIABLES= "X pixel", "Y pixel", "U pixel", "V pixel", "CHC" ZONE T="Pixel, Height=1004, Width=992" I=%d, J=%d, F=POINT\n',...
        fullfile(cd,sprintf('%s%04d.vec',fname,i)),r,c);
    fprintf(fid,'%s',header);
    fprintf(fid,'%8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n',cat(1,x,y,u(i,:),v(i,:),ones(1,length(x))));
    fclose(fid);
end
    