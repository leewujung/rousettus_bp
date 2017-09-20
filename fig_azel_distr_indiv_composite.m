% 2016 08 09  Plot distribution of azimuth and elevation for best-fitting
%             ellipse, for individual clicks and composite clicks

clear
warning off

if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/Dropbox/0_CODE/MATLAB/swtest');

    data_base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    usrn = getenv('username');
    if strcmp(usrn,'Wu-Jung')
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
        addpath('F:\Dropbox\0_CODE\MATLAB\swtest');
    else
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\swtest']);
    end

    if strcmp(usrn,'Wu-Jung')
        data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    else
        data_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
        save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
    end
end

results_path = 'analysis_results_figs';

indiv_data_path = 'rotate_all_click_20160721';
indiv_data_file = dir(fullfile(data_base_path,results_path,indiv_data_path,'*.mat'));

comp_data_path = 'multifreq_composite_click_fit_elps_20160808_batall_bin10_th2';
comp_data_file = dir(fullfile(data_base_path,results_path,comp_data_path,'*.mat'));
comp_data_file = comp_data_file.name;

ss = strsplit(comp_data_path,'_');
threshold = str2double(ss{end}(3:end));  % only use averaged results if number of measurement > threshold
binsize = str2double(ss{end-1}(4:end));  % bin size in azimuth and elevation
bat_num = ss{end-2}(4:end);

% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

title_text = sprintf('EL/AZ, indiv, comp: bat=%s, bin=%d, th=%d',...
    bat_num,binsize,threshold);
save_fname = sprintf('%s_indiv_comp_bat%s_bin%d_th%d',...
    script_name,bat_num,binsize,threshold);


for iF=1:length(indiv_data_file)
    B_data = load(fullfile(data_base_path,results_path,...
            indiv_data_path,indiv_data_file(iF).name));
    [ar_tmp,elps_xy_tmp] = get_ellipse_ar(B_data.rot_elpctr_tilt.E);
    ar(iF,1) = ar_tmp;
    elps_x_max(iF,1) = range(elps_xy_tmp(:,1))/2;
    elps_y_max(iF,1) = range(elps_xy_tmp(:,2))/2;
    e(iF,:) = B_data.rot_elpctr_tilt.E.e;
    a0(iF,:) = B_data.rot_elpctr_tilt.E.a0;
    b0(iF,:) = B_data.rot_elpctr_tilt.E.b0;
    theta(iF,:) = B_data.rot_elpctr_tilt.E.theta;
end

B_comp = load(fullfile(data_base_path,results_path,comp_data_path,comp_data_file));
freq_wanted = 35e3;
[~,freq_idx] = min(abs(B_comp.param.freq-freq_wanted));
[~,comp_elps_xy_tmp] = get_ellipse_ar(B_comp.left(freq_idx).rot_elpctr_tilt.E);
comp_left_elps_x_max = range(comp_elps_xy_tmp(:,1))/2;
comp_left_elps_y_max = range(comp_elps_xy_tmp(:,2))/2;
[~,comp_elps_xy_tmp] = get_ellipse_ar(B_comp.right(freq_idx).rot_elpctr_tilt.E);
comp_right_elps_x_max = range(comp_elps_xy_tmp(:,1))/2;
comp_right_elps_y_max = range(comp_elps_xy_tmp(:,2))/2;

% Convert from x-y to az-el angles
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

% Revmoe points out of xy map projection limit
[xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
xxx = xxx(1:2);
yyy = yyy(3:4);
xy_lim = [xxx(1) xxx(2) yyy(1) yyy(2)];
idx_bad = elps_x_max<xy_lim(1) | elps_x_max>xy_lim(2) |...
               elps_y_max<xy_lim(3) | elps_y_max>xy_lim(4);

% Inverse map projection
[el_data,az_data] = minvtran(mstruct,...
    elps_x_max(~idx_bad),elps_y_max(~idx_bad));
[el_comp_left,az_comp_left] = minvtran(mstruct,...
    comp_left_elps_x_max,comp_left_elps_y_max);
[el_comp_right,az_comp_right] = minvtran(mstruct,...
    comp_right_elps_x_max,comp_right_elps_y_max);

% Plot AZ/EL scatter plot
fig_az_el = figure;
plot(az_data(~idx_bad),el_data(~idx_bad),'.');
hold on
plot(az_comp_left,el_comp_left,'^');
plot(az_comp_right,el_comp_right,'s');
xlabel('Azimuth');
ylabel('Elevation');
legend('Individuals','Composite left','Composite right');
grid
axis equal
axis([0 90 0 90])
title(title_text);

saveas(fig_az_el,...
    fullfile(save_path,[save_fname,'_scatter.fig']),'fig');
saveSameSize(fig_az_el,'file',fullfile(save_path,[save_fname,'_scatter.png']),...
    'format','png','renderer','painters');

% Plot aspect ratio histogram
fig_ar = figure;
hd = histogram(el_data(~idx_bad)./az_data(~idx_bad),...
    0:0.1:4,'normalization','probability');
hold on
hleft = plot([1 1]*el_comp_left/az_comp_left,[0 1]);
hright = plot([1 1]*el_comp_right/az_comp_right,[0 1]);
xlabel('Ratio el/az');
ylabel('Relative frequency');
ylim([0 0.18])
title(title_text)

saveas(fig_ar,...
    fullfile(save_path,[save_fname,'_ar.fig']),'fig');
saveSameSize(fig_ar,'file',fullfile(save_path,[save_fname,'_ar.png']),...
    'format','png','renderer','painters');
