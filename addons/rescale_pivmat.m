function data = rescale_pivmat(data,m2p,dt)
% RESCALE_PIVMAT(PIVMAT,METER-TO-PIXEL,DT) rescales the PIV data in CIL MAT format
% Basically it does the following set:
% PIVMAT.X = PIVMAT.X * METER-TO-PIXEL
% PIVMAT.U = PIVMAT.U * METER-TO-PIXEL/DT
%
% Example:
% re1500 = load('/Users/alex/Desktop/re1500.mat')
% re1500 = rescale_pivmat(re1500,0.08/1888, 30e-3); % Cavity data: 0.08 m = 1888
% pixels and 30 milliseconds is the DT
% showf(averf(cilmat2pivmat(re1500)))
%


% Author: Alex Liberzon (alex dot liberzon at gmail dot com)
% Copyright (c) 1998-2012 OpenPIV group
% See the file license.txt for copying permission.



for i = 1:length(data)
    
    data(i).x = (data(i).x - min(data(i).x)) * m2p;
    data(i).y = (data(i).y - min(data(i).y)) * m2p;
    data(i).vx = data(i).vx * m2p/dt;
    data(i).vy = data(i).vy * m2p/dt;
    data(i).unitx = 'm';
    data(i).unity = 'm';
    data(i).unitvx = 'm/s';
    data(i).unitvy = 'm/s';
    
end