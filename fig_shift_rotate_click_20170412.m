% 2016 08 15  Plotting the procedure to shift/rotate for aligning clicks
% 2017 03 01  Replot using new results with corrected bat_to_mic_angle

clear

if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');
    
    data_base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\mtit');
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end

% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));

results_path = 'analysis_results_figs';
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Set params
bat_proc_path = 'rotate_all_click_2010308';
bat_num = '36134';
trial_num = 2;
click_num = 19;
freq_wanted = 35*1e3;
interp_opt = 'rb_rbf';
% cvec = 0:-3:-30;
cvec = [0:-3:-30];

% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

% Load rotated data
bat_proc_file = sprintf('%s_%s_%02d_c%02d_rotated.mat',...
    bat_proc_path,bat_num,trial_num,click_num);
save_fname = sprintf('%s_%s_%02d_c%02d_%2dkHz',...
    script_name,bat_num,trial_num,click_num,freq_wanted);

V = load(fullfile(data_base_path,results_path,bat_proc_path,bat_proc_file));
suptitle_text = sprintf('Bat %s, Trial %d, Call #%02d',...
    bat_num,trial_num,click_num);

% Movement of [az,el] of each mic
fig_azel = plot_indiv_click_azel_movement(V);
title(suptitle_text,'fontsize',14)
saveSameSize(fig_azel,'file',fullfile(save_path,[save_fname,'_azel.png']),...
    'format','png','renderer','painters');
close(fig_azel)

% Plot the shift/tilt procedure
V.raw.vq_norm(V.raw.vq_norm<-30) = -30;
V.rot_max.vq_norm(V.rot_max.vq_norm<-30) = -30;
V.rot_elpctr.vq_norm(V.rot_elpctr.vq_norm<-30) = -30;
V.rot_elpctr_tilt.vq_norm(V.rot_elpctr_tilt.vq_norm<-30) = -30;
fig_fit = plot_indiv_click_rotate(V,cvec,mstruct);
mtit(fig_fit,suptitle_text,'fontsize',14)
saveas(fig_fit,fullfile(save_path,[save_fname,'_rot.fig']));
saveSameSize(fig_fit,'file',fullfile(save_path,[save_fname,'_rot.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_rot.eps']));
close(fig_fit)




