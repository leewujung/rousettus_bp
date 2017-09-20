% 2015 12 28  Plots for multi-freq composite clicks
% 2016 04 20  Separate this part out from multifreq_composite_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 07 23  Plot for paper; best-fitting ellipse for composite click
% 2016 08 08  Plot for paper; revised to use newly calculated results from
%             multifreq_composite_click_20160808.m
% 2016 10 25  Update for version 1025 with out-of-bound points

clear

if isunix
else
    usrn = getenv('username');
    if strcmp(usrn,'Wu-Jung')
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
        addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
    else
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
    end
    
    if strcmp(usrn,'Wu-Jung')
        data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    else
        data_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
        save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    end
end

% Set path and params
results_path = 'analysis_results_figs';
data_path = 'multifreq_composite_click_20161025_batall_bin10_th0';
data_file = sprintf('%s_merged_clicks.mat',data_path);

ss = strsplit(data_path,'_');
threshold = str2double(ss{end}(3:end));  % only use averaged results if number of measurement > threshold
binsize = str2double(ss{end-1}(4:end));  % bin size in azimuth and elevation
bat_num = ss{end-2}(4:end);

it_shift_th = 0.005;
azel_bnd = [300,150];
cgrey = 200*ones(1,3)/255;

% Save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_bat%s_bin%d_th%d',script_name,bat_num,binsize,threshold);
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end
save_mat_fname = sprintf('%s.mat',script_name);

% Load data
D = load(fullfile(data_base_path,results_path,data_path,data_file));

% Save params
S.param.data_path = data_path;
S.param.data_file = data_file;
S.param.it_shift_th = it_shift_th;
S.param.azel_bnd = azel_bnd;
S.param.freq = D.param.freq_wanted;

num_freq = length(D.param.freq_wanted);
for iF=1:num_freq
    freq_wanted = D.param.freq_wanted(iF);
    
    % params for figure
    vq_norm_min = -27;
    cvec = 0:-3:(floor(vq_norm_min/3)-1)*3;
    cvec_min_idx = find(cvec-vq_norm_min<0,1,'first');
    map_proj = D.param.map.map_projection;
    
    S.cvec = cvec;
    S.param.map = D.param.map;
    
    % Right click =============
    title_text = sprintf('Averaged right clicks, bat %s, %d kHz, th=%d, bin=%ddeg',...
        bat_num,freq_wanted/1e3,D.param.composite_threshold,D.param.composite_binsize);
    save_fname = sprintf('%s_%dkHz_right',script_name,freq_wanted/1e3);
    
    azq = D.averaged_composite.right.interp(iF).azq_avg;
    elq = D.averaged_composite.right.interp(iF).elq_avg;
    vq_norm = D.averaged_composite.right.interp(iF).vq_norm_avg;
    
    [raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(azq,elq,vq_norm,map_proj,it_shift_th,azel_bnd);
    S.right(iF).raw = raw;
    S.right(iF).rot_max = rot_max;
    S.right(iF).rot_elpctr = rot_elpctr;
    S.right(iF).rot_elpctr_tilt = rot_elpctr_tilt;
    
    % Plot rotation procedure
    fig_fit = plot_indiv_click_rotate(S.right(iF),cvec,D.param.map.mstruct);
    suptitle(title_text);
    saveSameSize(fig_fit,'file',...
        fullfile(save_path,save_fname),...
        'format','png','renderer','painters');
    close(fig_fit)
    
    % Left click =============
    title_text = sprintf('Averaged left clicks, bat %s, %d kHz, th=%d, bin=%ddeg',...
        bat_num,freq_wanted/1e3,D.param.composite_threshold,D.param.composite_binsize);
    save_fname = sprintf('%s_%dkHz_left',script_name,freq_wanted/1e3);
    
    azq = D.averaged_composite.left.interp(iF).azq_avg;
    elq = D.averaged_composite.left.interp(iF).elq_avg;
    vq_norm = D.averaged_composite.left.interp(iF).vq_norm_avg;
    
    [raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(azq,elq,vq_norm,map_proj,it_shift_th,azel_bnd);
    S.left(iF).raw = raw;
    S.left(iF).rot_max = rot_max;
    S.left(iF).rot_elpctr = rot_elpctr;
    S.left(iF).rot_elpctr_tilt = rot_elpctr_tilt;
    
    % Plot rotation procedure
    fig_fit = plot_indiv_click_rotate(S.left(iF),cvec,D.param.map.mstruct);
    suptitle(title_text);
    saveSameSize(fig_fit,'file',...
        fullfile(save_path,save_fname),...
        'format','png','renderer','painters');
    close(fig_fit)
end

save(fullfile(save_path,[script_name,'.mat']),'-struct','S');


% Plot best-fitting ellipse param
for iF=1:num_freq
    [ar_left(iF),xy] = get_ellipse_ar(S.left(iF).rot_elpctr_tilt.E);
    x_max = range(xy(:,1))/2;
    y_max = range(xy(:,2))/2;
    [el,az] = minvtran(S.param.map.mstruct,x_max,y_max);
    azel_ar_left(iF) = el/az;
    
    [ar_right(iF),xy] = get_ellipse_ar(S.right(iF).rot_elpctr_tilt.E);
    x_max = range(xy(:,1))/2;
    y_max = range(xy(:,2))/2;
    [el,az] = minvtran(S.param.map.mstruct,x_max,y_max);
    azel_ar_right(iF) = el/az;
end


figure
plot(S.param.freq/1e3,azel_ar_left,'-o');
hold on
plot(S.param.freq/1e3,azel_ar_right,'-o');
legend('Left','Right')
grid
xlabel('Frequency (kHz)');
ylabel('EL/AZ aspect ratio');
title_text = sprintf('EL/AZ, bat %s, th=%d, bin=%ddeg',...
        bat_num,D.param.composite_threshold,D.param.composite_binsize);
title(title_text);
saveSameSize(gcf,'file',...
    fullfile(save_path,sprintf('%s_azel_ratio.png',script_name)),...
    'format','png','renderer','painters');
close




