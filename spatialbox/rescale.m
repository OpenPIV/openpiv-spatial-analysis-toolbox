function data = rescale(data,m2p,dt)
% RESCALE(CILMAT,METER-TO-PIXEL,DT) rescales the PIV data in CIL MAT format
% Basically it does the following set:
% CILMAT.X = CILMAT.X * METER-TO-PIXEL
% CILMAT.U = CILMAT.U * METER-TO-PIXEL/DT
%
% Example:
% re1500 = load('/Users/alex/Desktop/re1500.mat')
% re1500 = rescale(re1500,0.08/1888, 30e-3); % Cavity data: 0.08 m = 1888
% pixels and 30 milliseconds is the DT
% showf(averf(cilmat2pivmat(re1500)))
% 

data.x = data.x * m2p;
data.y = data.y * m2p;
data.u = data.u * m2p/dt;
data.v = data.v * m2p/dt;
data.uf = data.uf * m2p/dt;
data.vf = data.vf * m2p/dt;

data.xUnits = 'm'
data.velUnits = 'm/s'