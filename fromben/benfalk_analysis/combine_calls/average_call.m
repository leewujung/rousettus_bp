function [interp_avg,bin_avg] = average_call(az,el,call_dB,binsize,mproj,threshold)
% Merge data from all calls, only using bins containing >1 measurements
%
% INPUT
%   az        azimuth angle [deg]
%   el        elevation angle [deg]
%   call_dB   call measurement [dB]
%   binsize   size of angle bin [deg]
%   mproj     map projection to use
%   threshold   threshold to average or discard measurement
%               (only use averaged values when then number of measurements > threshold)
%
% OUTPUT
%   bin_avg      structure containing *binned* averaged output with azimuth,
%                elevation, and call_dB
%   interp_avg   structure contatining *interpolated* results with azimuth,
%                elevation, call_dB, and map projected x-y coord
%   all output angles are in [deg]
%
% Wu-Jung Lee
% 2015 12 03

[n,az_edge,el_edge] = histcounts2(az,el,-180:binsize:181,-90:binsize:91);
nidx = find(n>threshold);
[az_edge_idx,el_edge_idx] = ind2sub(size(n),nidx);
avg_call_dB = nan(size(n));
for iN = 1:length(nidx)
    avg_call_dB(nidx(iN)) = nanmean(call_dB( az>=az_edge(az_edge_idx(iN)) & az<az_edge(az_edge_idx(iN)+1) &...
                                             el>=el_edge(el_edge_idx(iN)) & el<el_edge(el_edge_idx(iN)+1)));
end

idx_avg = find(~isnan(avg_call_dB));
[elmesh,azmesh] = meshgrid(el_edge(1:end-1)+mean(diff(el_edge))/2,az_edge(1:end-1)+mean(diff(az_edge))/2);
az_avg = azmesh(idx_avg);
el_avg = elmesh(idx_avg);
avg_call_dB = avg_call_dB(idx_avg);
[~,vq_norm_avg,azq_avg,elq_avg] = interp_bp(az_avg/180*pi,el_avg/180*pi,avg_call_dB,'natural');

mstruct = defaultm(mproj);
mstruct = defaultm(mstruct);
[xq_avg,yq_avg] = mfwdtran(mstruct,elq_avg/pi*180,azq_avg/pi*180);

interp_avg.vq_norm_avg = vq_norm_avg;
interp_avg.azq_avg = azq_avg/pi*180;
interp_avg.elq_avg = elq_avg/pi*180;
interp_avg.xq_avg = xq_avg;
interp_avg.yq_avg = yq_avg;
bin_avg.avg_call_dB = avg_call_dB;
bin_avg.az_avg = az_avg;
bin_avg.el_avg = el_avg;



