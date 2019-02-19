function mat2vec(x,z,u,w,fname,varargin)
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
[r,c,k] = size(u);
x = reshape(permute(x,[2 1]),[1, r*c]);
z = reshape(permute(z(end:-1:1,:),[2 1]),[1, r*c]);
u = reshape(permute(u,[2 1 3]),[r*c k]);
w = reshape(permute(w,[2 1 3]),[r*c k]);
% x = reshape(permute(x(end:-1:1,:),[2 1]),[1, r*c]);
% z = reshape(permute(z(end:-1:1,:),[2 1]),[1, r*c]);
% u = reshape(permute(u(end:-1:1,:,:),[2 1 3]),[r*c k]);
% w = reshape(permute(w(end:-1:1,:,:),[2 1 3]),[r*c k]);



for i = 1:k
    fid = fopen(sprintf('%s%04d.vec',fname,i),'w');
    header = sprintf('TITLE="%s" VARIABLES= "X pixel", "Y pixel", "U pixel", "V pixel", "CHC" ZONE T="Pixel, Height=1004, Width=992" I=%d, J=%d, F=POINT\n',...
        fullfile(cd,sprintf('%s%04d.vec',fname,i)),r,c);
    fprintf(fid,'%s',header);
    fprintf(fid,'%8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n',cat(1,x,z,u(:,i)',w(:,i)',ones(1,length(x))));
    fclose(fid);
end
    