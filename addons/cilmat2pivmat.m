function [varargout] = cilmat2pivmat(varargin)
% CILMAT2PIVMAT converts the CIL Spatialbox MAT files into PIVMAT
% compatible MAT files
% main differences: CILMAT has velocity in 3D matrices, and PIVMAT in
% large structure, e.g.: 
% u = struct 1 x 1, in which matrix dudx is R x C x N (fields)
% v = struct 1 x N (fields)
% R = size in rows
% C = size in cols
% main translation is the loop:
% for i = 1:size(u.uf,3) % for all fields
%   newstruct(i).uf = u.uf(:,:,i); 
% end
% 

% Author: Alex Liberzon (alex dot liberzon at gmail dot com)
% Copyright (c) 1998-2012 OpenPIV group
% See the file license.txt for copying permission.


u = varargin{1};

v = struct('x',[],'y',[],'vx',[],'vy',[],'unitx',[],'unity',[],...
    'unitvx',[],'unitvy',[],'source',[],'name',[],...
    'namex','x','namey','y','namevx','u_x','namevy','u_y',...
    'history',{{'cilmat2pivmat'}},'ysign','Y axis upward','setname','test');
v = repmat(v,[size(u.uf,3),1]);
for i = 1:size(u.uf,3)
    v(i).x = u.x(1,:)
    v(i).y = u.y(:,1)'; 
    v(i).vx = u.u(:,:,i).';
    v(i).vy = u.v(:,:,i).';
%     v(i).choice = [];
    v(i).unitx = u.xUnits;
    v(i).unity = u.xUnits; 
    v(i).unitvx = u.velUnits;
    v(i).unitvy = u.velUnits;
    v(i).source = 'Insight,CIL';
    v(i).name = u.path;
end
    
varargout{1} = v;

