% 2016 08 09  Plot distribution of azimuth and elevation for best-fitting
%             ellipse, for individual clicks and composite clicks
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 04 12  Replot using new results with corrected bat_to_mic_angle
% 2017 09 23  Also plot the azel distribution of Monte Carlo simulated model
%             beampattern


clear

if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/swtest');
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\swtest');
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end


% Path to data and model
results_path = 'analysis_results_figs';

indiv_data_path = 'rotate_all_click_2010308';
indiv_data_file = dir(fullfile(data_base_path,results_path,indiv_data_path,'*.mat'));

comp_data_path = 'multifreq_composite_click_fit_elps_20170412_batall_bin10_th0';
comp_data_file = dir(fullfile(data_base_path,results_path,comp_data_path,'*.mat'));
comp_data_file = comp_data_file.name;

piston_path = 'model_piston_proj_20170412_azel_std1.0';
piston_file = dir(fullfile(data_base_path,results_path,piston_path,'*.mat'));

phase_path = 'model_bp_proj_RaColony_multifreq_20170921_std1.0';
phase_file = dir(fullfile(data_base_path,results_path,phase_path,'*.mat'));


% File name stuff
ss = strsplit(comp_data_path,'_');
threshold = str2double(ss{end}(3:end));  % only use averaged results if number of measurement > threshold
binsize = str2double(ss{end-1}(4:end));  % bin size in azimuth and elevation
bat_num = ss{end-2}(4:end);

noise_std = piston_path(end-2:end);

param.noise_std = noise_std;
param.binsize = binsize;
param.threshold = threshold;
param.bat_num = bat_num;

% Other params
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

[xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
xxx = xxx(1:2);
yyy = yyy(3:4);
xy_lim = [xxx(1) xxx(2) yyy(1) yyy(2)];

freq_wanted = 35e3;

param.freq_wanted = freq_wanted;
param.map.map_proj = map_proj;
param.map.mstruct = mstruct;

% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

title_text = sprintf('EL/AZ, indiv, comp: bat=%s, bin=%d, th=%d, nstd=%s',...
    bat_num,binsize,threshold,noise_std);
save_fname = sprintf('%s_bat%s_bin%d_th%d_nstd%s',...
    script_name,bat_num,binsize,threshold,noise_std);


% Get az/el params for all data and models
for iF=1:length(indiv_data_file)

    % individual click data
    B_data = load(fullfile(data_base_path,results_path,...
                           indiv_data_path,indiv_data_file(iF).name));
    data(iF) = get_azel_params(B_data);

    % piston model
    B_piston = load(fullfile(data_base_path,results_path,...
                             piston_path,piston_file(iF).name));
    piston(iF) = get_azel_params(B_piston);

    % phased array model
    B_phase = load(fullfile(data_base_path,results_path,...
                             phase_path,phase_file(iF).name));
    phase(iF) = get_azel_params(B_phase);
end

% Composite measured clicks
B_comp = load(fullfile(data_base_path,results_path,comp_data_path,comp_data_file));
[~,freq_idx] = min(abs(B_comp.param.freq-freq_wanted));
comp_left = get_azel_params(B_comp.left(freq_idx));
comp_right = get_azel_params(B_comp.right(freq_idx));


% Remove points out of xy map projection limit
data_idx_bad = [data(:).elps_x]<xy_lim(1) | [data(:).elps_x]>xy_lim(2) |...
               [data(:).elps_y]<xy_lim(3) | [data(:).elps_y]>xy_lim(4);
piston_idx_bad = [piston(:).elps_x]<xy_lim(1) | [piston(:).elps_x]>xy_lim(2) |...
                 [piston(:).elps_y]<xy_lim(3) | [piston(:).elps_y]>xy_lim(4);
phase_idx_bad = [phase(:).elps_x]<xy_lim(1) | [phase(:).elps_x]>xy_lim(2) |...
                [phase(:).elps_y]<xy_lim(3) | [phase(:).elps_y]>xy_lim(4);


% Convert from x-y to az-el angles -- inverse map projection
[el_data,az_data] = minvtran(mstruct,[data(~data_idx_bad).elps_x]',...
                                     [data(~data_idx_bad).elps_y]');
[el_piston,az_piston] = minvtran(mstruct,[piston(~piston_idx_bad).elps_x]',...
                                         [piston(~piston_idx_bad).elps_y]');
[el_phase,az_phase] = minvtran(mstruct,[phase(~phase_idx_bad).elps_x]',...
                                       [phase(~phase_idx_bad).elps_y]');
[el_comp_left,az_comp_left] = minvtran(mstruct,comp_left.elps_x,comp_left.elps_y);
[el_comp_right,az_comp_right] = minvtran(mstruct,comp_right.elps_x,comp_right.elps_y);



% Statistical tests
% --- KS 2-sample test
[ks.h_data_phase,ks.p_data_phase] = ...
    kstest2(el_data(~data_idx_bad)./az_data(~data_idx_bad),...
            el_phase(~phase_idx_bad)./az_phase(~phase_idx_bad));
[ks.h_data_piston,ks.p_data_piston] = ...
    kstest2(el_data(~data_idx_bad)./az_data(~data_idx_bad),...
            el_piston(~piston_idx_bad)./az_piston(~piston_idx_bad));
% Mann-Whitney U-test
[mw.h_data_phase,mw.p_data_phase,mw.stat_data_phase] = ...
    ranksum(el_data(~data_idx_bad)./az_data(~data_idx_bad),...
            el_phase(~phase_idx_bad)./az_phase(~phase_idx_bad));
[mw.h_data_piston,mw.p_data_piston,mw.stat_data_phase] = ...
    ranksum(el_data(~data_idx_bad)./az_data(~data_idx_bad),...
            el_piston(~piston_idx_bad)./az_piston(~piston_idx_bad));
stat.ks = ks;
stat.mw = mw;



% Save results
save(fullfile(save_path,[save_fname,'_results.mat']),...
     'param','data*','piston*','phase*','el_*','az_*','stat');


% Plot AZ/EL scatter plot
fig_az_el = figure;
corder = get(gca,'colororder');
hm = plot(az_piston,el_piston,'ko','markersize',2);
hold on
hp = plot(az_phase,el_phase,'.','color',corder(3,:),'markersize',7);
hd = plot(az_data,el_data,'.','color',corder(1,:),'markersize',7);
hcl = plot(az_comp_left,el_comp_left,'^','color',corder(2,:),'linewidth',2);
hcr = plot(az_comp_right,el_comp_right,'s','color',corder(2,:),'linewidth',2);
xlabel('Azimuth');
ylabel('Elevation');
legend([hd,hm,hp,hcl,hcr],'Indiv clicks','Piston','Phased array',...
       'Composite left','Composite right');
grid
axis equal
axis([0 90 0 90])
title(title_text);

saveas(fig_az_el,...
    fullfile(save_path,[save_fname,'_scatter.fig']),'fig');
saveSameSize_res(fig_az_el,150,'file',fullfile(save_path,[save_fname,'_scatter.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_scatter.eps']));



% Plot aspect ratio histogram
fig_ar = figure('position',[50,50,350,450]);
subplot(311)
hd = histogram(el_data./az_data,0:0.1:4,'normalization','probability');
hold on
hleft = plot([1 1]*el_comp_left/az_comp_left,[0 1]);
hright = plot([1 1]*el_comp_right/az_comp_right,[0 1]);
legend([hd,hleft,hright],{'Indiv clicks','Composite left','Composite right'});
ylim([0 0.18])
title(title_text)
subplot(312)
hm = histogram(el_piston./az_piston,...
               0:0.1:4,'normalization','probability');
legend('Piston')
subplot(313)
hp = histogram(el_phase./az_phase,...
               0:0.1:4,'normalization','probability');
legend('Phased array')
xlabel('Ratio el/az');
ylabel('Relative frequency');
ylim([0 0.35])

saveas(fig_ar,...
    fullfile(save_path,[save_fname,'_ar.fig']),'fig');
saveSameSize_res(fig_ar,150,'file',fullfile(save_path,[save_fname,'_ar.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_ar.eps']));



