% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 08 09  Plot for paper;
%             revised from multifreq_composite_click_cntr_20160808
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 03 08  Re-processed all files due to errors in mic_bp beampattern
%             compensation and degF to degC conversion.
% 2017 09 23  Plot beam center along with multifreq contour
%             The beam center finding routines are from `fig_composite_click_avg_bp_20170920`
% 2017 09 23  Use the beam center finding routines from
%             `fig_composite_click_avg_bp_20170920` to calculate and add beam
%             center info to the composite bp output from `multifreq_composite_click_20170308`

clear

if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');

    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end


% Set up various paths
results_path = 'analysis_results_figs';
data_path = 'multifreq_composite_click_20170308_batall_bin10_th0';

% Load data
data_file = sprintf('%s_merged_clicks.mat',data_path);
D = load(fullfile(data_base_path,results_path,data_path,data_file));

num_freq = length(D.param.freq_wanted);
map_proj = 'eckert4';
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
D.map.map_projection = map_proj;
D.map.mstruct = mstruct;

for iF=1:num_freq
    
    % Find beam center based on energy
    % --- left composite click
    xx = D.averaged_composite.left.interp(iF).vq_norm_avg(:);
    [left_max_val,left_max_idx] = max(xx);
    D.bpctr{iF}.left.max.el = D.averaged_composite.left.interp(iF).elq_avg(left_max_idx);
    D.bpctr{iF}.left.max.az = D.averaged_composite.left.interp(iF).azq_avg(left_max_idx);
    xx(isnan(xx)) = -inf;
    [~,left_sort_idx] = sort(xx,'descend');
    ii = xx(left_sort_idx)>-1;
    D.bpctr{iF}.left.top.el = mean(D.averaged_composite.left.interp(iF).elq_avg(left_sort_idx(ii)));
    D.bpctr{iF}.left.top.az = mean(D.averaged_composite.left.interp(iF).azq_avg(left_sort_idx(ii)));
    % --- right composite click
    xx = D.averaged_composite.right.interp(iF).vq_norm_avg(:);
    [right_max_val,right_max_idx] = max(xx);
    D.bpctr{iF}.right.max.el = D.averaged_composite.right.interp(iF).elq_avg(right_max_idx);
    D.bpctr{iF}.right.max.az = D.averaged_composite.right.interp(iF).azq_avg(right_max_idx);
    xx(isnan(xx)) = -inf;
    [~,right_sort_idx] = sort(xx,'descend');
    ii = xx(right_sort_idx)>-1;
    D.bpctr{iF}.right.top.el = mean(D.averaged_composite.right.interp(iF).elq_avg(right_sort_idx(ii)));
    D.bpctr{iF}.right.top.az = mean(D.averaged_composite.right.interp(iF).azq_avg(right_sort_idx(ii)));

    % Fit ellipse
    % --- left composite click
    [left_raw,left_rot_max,left_rot_elpctr,left_rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(D.averaged_composite.left.bin(iF).az_avg,...
                                  D.averaged_composite.left.bin(iF).el_avg,...
                                  D.averaged_composite.left.bin(iF).avg_call_dB,...
                                  map_proj,0.005);
    [left_el_ectr,left_az_ectr] = minvtran(mstruct,left_rot_max.E.x0,left_rot_max.E.y0);  % inverse map projection
    [D.bpctr{iF}.left.ectr.el,D.bpctr{iF}.left.ectr.az] = rotatem(left_el_ectr,left_az_ectr,...
                                    [D.bpctr{iF}.left.max.el,D.bpctr{iF}.left.max.az],...
                                    'inverse','degrees');
    % --- right composite click
    [right_raw,right_rot_max,right_rot_elpctr,right_rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(D.averaged_composite.right.bin(iF).az_avg,...
                                  D.averaged_composite.right.bin(iF).el_avg,...
                                  D.averaged_composite.right.bin(iF).avg_call_dB,...
                                  map_proj,0.005);
    [right_el_ectr,right_az_ectr] = minvtran(mstruct,right_rot_max.E.x0,right_rot_max.E.y0);  % inverse map projection
    [D.bpctr{iF}.right.ectr.el,D.bpctr{iF}.right.ectr.az] = rotatem(right_el_ectr,right_az_ectr,...
                                    [D.bpctr{iF}.right.max.el,D.bpctr{iF}.right.max.az],...
                                    'inverse','degrees');

end

save(fullfile(data_base_path,results_path,data_path,data_file),'-struct','D');



