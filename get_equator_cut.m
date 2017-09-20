function [vq_xy,x,y] = get_equator_cut(dataM,call_dB,mstruct)
% Get interpolated value along equator (elevation=0)

% INPUT
%    M         structure that contains x,y,xq,yq
%    call_dB   measurements at (x,y)

% Wu-Jung Lee | leewujung@gmail.com
% 2015 12 24

% Get interpolated values along equator
az_cut = -180:3:180;
[x,~] = mfwdtran(mstruct,zeros(1,length(az_cut)),az_cut);
[x,y] = meshgrid(x,0);
vq_xy = rbfinterp([x;y],rbfcreate([dataM.x(:)';dataM.y(:)'],call_dB-max(call_dB),'RBFFunction','multiquadrics'));
vq_xy = reshape(vq_xy,size(x));

% Set values outside of boundary to NaN
k = boundary(dataM.x(:),dataM.y(:),0);  % outer boundary of all measured points
[in,on] = inpolygon(x,y,dataM.x(k),dataM.y(k));
in = in|on;
vq_xy(~in) = NaN;
