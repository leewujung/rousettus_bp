function [vq,vq_norm,azq,elq] = interp_bp(az,el,mic_dB,method)
% az,el    azimuth and elevation location for each mic [radian]
% mic_dB   corresponding mic recording amplitude [dB]
% method   interpolation method, 'natural' or 'rbf'

% Wu-Jung Lee | leewujung@gmail.com
% 2016/10/25  Due to changes in shift_rotate_bp, now interp_bp may have
%             az/el = NaN

% Get data
az = az(:);
el = el(:);
mic_dB = mic_dB(:);

% Below added on 2016/10/25
idx_nan = isnan(az) | isnan(el) | isnan(mic_dB);
az(idx_nan) = [];
el(idx_nan) = [];
mic_dB(idx_nan) = [];

% Interpolation
maxref = max(mic_dB);
[azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
if strcmp(method,'natural')  % natural neighbor interpolation
    vq = griddata(az,el,mic_dB,azq,elq,'natural');
elseif strcmp(method,'rbf')  % radial basis function interpolation
    vq = rbfinterp([azq(:)';elq(:)'],rbfcreate([az';el'],mic_dB','RBFFunction','multiquadrics'));
    vq = reshape(vq,size(azq));
end
vq_norm = vq-maxref;

% Find boundary
idx_good = ~isnan(az) & ~isnan(el);
az = az(idx_good);
el = el(idx_good);

k = boundary(az,el,0);  % outer boundary of all measured points
[in,on] = inpolygon(azq,elq,az(k),el(k));
in = in|on;

% Set values outside of az-el boundary to NaN
vq(~in) = NaN;
vq_norm(~in) = NaN;

