% 2016 08 08  Assemble composite clicks from simulated mic receptions
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 09 23  Enable assembling composite click from multi-freq simulation

clear

usrn = getenv('username');
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    model_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling/';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    model_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
end

results_path = 'analysis_results_figs';
simu_data_path = 'model_bp_proj_RaColony_multifreq_20170921_std1.0';
ss = strsplit(simu_data_path,'_');
model_shape = ss{4};
model_calc_date = ss{6};
noise_std = ss{end};
simu_file = dir(fullfile(data_base_path,results_path,simu_data_path,'*.mat'));

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

A.param.data_base_path = data_base_path;
A.param.simu_data_path = simu_data_path;
A.param.simu_data_file = simu_file;

% Composite click params
threshold = 0;  % only use averaged results if number of measurement > threshold
binsize = 10;  % bin size in azimuth and elevation

A.param.composite_threshold = threshold;
A.param.composite_binsize = binsize;

% Load 1 file to get mtx sizes
D = load(fullfile(data_base_path,results_path,simu_data_path,simu_file(1).name));  % rotated data
num_ch = size(D.v_mic,1);
num_freq = length(D.freq.all);

A.param.num_ch = num_ch;
A.param.map = D.map;

% Set save_fname
save_fname = sprintf('%s_%s_bin%d_th%d',...
    script_name,noise_std,binsize,threshold);

% Merge all data
click_side_all = nan(length(simu_file),num_ch);
az_elpctr_tilt_all = nan(length(simu_file),num_ch);
el_elpctr_tilt_all = nan(length(simu_file),num_ch);
call_dB_all = nan(length(simu_file),num_ch,length(D.freq.all));
for iS=1:length(simu_file)

    D = load(fullfile(data_base_path,results_path,simu_data_path,simu_file(iS).name));  % rotated data

    % Merge all data
    if isempty(D.rot_elpctr_tilt)
        click_side_all(iS,:) = nan(1,num_ch);
        az_elpctr_tilt_all(iS,:) = nan(1,num_ch);
        el_elpctr_tilt_all(iS,:) = nan(1,num_ch);
        call_dB_all(iS,:,:) = nan(num_ch,num_freq);
    else
        click_side_all(iS,:) = D.click_side*ones(1,length(num_ch));
        az_elpctr_tilt_all(iS,:) = D.rot_elpctr_tilt.az;
        el_elpctr_tilt_all(iS,:) = D.rot_elpctr_tilt.el;
        call_dB_all(iS,:,:) = D.v_mic;
    end

end
A.click_side_all = click_side_all;
A.az_elpctr_tilt_all = az_elpctr_tilt_all;
A.el_elpctr_tilt_all = el_elpctr_tilt_all;
A.call_dB_all = call_dB_all;

% All data points for a single freq
click_side_merge = reshape(click_side_all,[],1);
az_elpctr_tilt_merge = reshape(az_elpctr_tilt_all,[],1);
el_elpctr_tilt_merge = reshape(el_elpctr_tilt_all,[],1);
call_dB_merge = reshape(call_dB_all,[],4);

idx_right_good = find(click_side_merge==1 & ~isnan(call_dB_merge(:,1)));
call_dB_right = call_dB_merge(idx_right_good,:);
az_right = az_elpctr_tilt_all(idx_right_good);
el_right = el_elpctr_tilt_all(idx_right_good);

idx_left_good = find(click_side_merge==0 & ~isnan(call_dB_merge(:,1)));
call_dB_left = call_dB_merge(idx_left_good,:);
az_left = az_elpctr_tilt_all(idx_left_good);
el_left = el_elpctr_tilt_all(idx_left_good);



% Interpolation using averaged data
for iF=1:num_freq
    fprintf('Merging freq %d kHz\n',D.freq.all(iF)/1e3);
    [interp_avg_right(iF),bin_avg_right(iF)] = ...
        average_call(az_right,el_right,call_dB_right(:,iF),binsize,'eckert4',threshold);
    [interp_avg_left(iF),bin_avg_left(iF)] = ...
        average_call(az_left,el_left,call_dB_left(:,iF),binsize,'eckert4',threshold);
end


% Save mat file
save(fullfile(save_path,[save_fname,'.mat']),'-struct','A');



% Plot params
suptitle_str_pre = sprintf('%s %s',model_shape,noise_std);
cgrey = 200*ones(1,3)/255;
cvec=-30:3:0;


% Plot all scatter samples at all freq
fig_scatter = figure('position',[200,50,800,950]);

for iFcount=1:num_freq
    if iFcount==1  % swap 25 and 35 kHz locations
        iF = 2;
    elseif iFcount==2
        iF = 1;
    else
        iF = iFcount;
    end

    fprintf('Plotting scatter at %d kHz\n',D.freq.all(iF)/1e3);

    subplot(num_freq,2,(iFcount-1)*2+1)  % left click
    axesm(D.map.map_projection);
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    scatterm(el_left,az_left,15,call_dB_left(:,iF),'filled')
    title(sprintf('%dkHz left',D.freq.all(iF)/1e3));
    tightmap
    axis off
    caxis(cvec([1 end]))
    %colorbar('location','eastoutside');

    subplot(num_freq,2,(iFcount-1)*2+2)  % right click
    axesm(D.map.map_projection);
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    scatterm(el_right,az_right,15,call_dB_right(:,iF),'filled')
    title(sprintf('%dkHz right',D.freq.all(iF)/1e3));
    tightmap
    axis off
    caxis(cvec([1 end]))
    %colorbar('location','eastoutside');
end

figure(fig_scatter)
saveSameSize(fig_scatter,'file',fullfile(save_path,[save_fname,'_scatter.png']),...
             'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_scatter.eps']));


% Plot reconstructed bp at all freq
fig_avg_bp_all = figure('position',[200,50,800,950]);

for iFcount=1:num_freq
    if iFcount==1  % swap 25 and 35 kHz locations
        iF = 2;
    elseif iFcount==2
        iF = 1;
    else
        iF = iFcount;
    end

    fprintf('Plotting scatter at %d kHz\n',D.freq.all(iF)/1e3);

    % left click
    plot_bp_simple(subplot(num_freq,2,(iFcount-1)*2+1),...
                   interp_avg_left(iF).azq_avg,interp_avg_left(iF).elq_avg,...
                   interp_avg_left(iF).vq_norm_avg,'eckert4');
    title('Reconstructed bp: left');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    % right click
    plot_bp_simple(subplot(num_freq,2,(iFcount-1)*2+2),...
                   interp_avg_right(iF).azq_avg,interp_avg_right(iF).elq_avg,...
                   interp_avg_right(iF).vq_norm_avg,'eckert4');
    title('Reconstructed bp: right');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');
end

figure(fig_avg_bp_all)
saveSameSize(fig_avg_bp_all,'file',fullfile(save_path,[save_fname,'_avg_bp_all.png']),...
             'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_avg_bp_all.eps']));



% Plot model bp and reconstructed bp for comparison
for iF = 1:num_freq

    fprintf('Plotting reconstructed bp at %d kHz',D.freq.all(iF)/1e3);

    fig_bp_cmp = figure('position',[200,60,1200,900]);

    % model bp: left
    plot_bp_simple(subplot(221),-D.BP(iF).az,D.BP(iF).el,D.BP(iF).pp_plot, ...
                   D.map.map_projection);
    title('model: left');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    % model bp: right
    plot_bp_simple(subplot(222),D.BP(iF).az,D.BP(iF).el,D.BP(iF).pp_plot, ...
                   D.map.map_projection);
    title('model: right');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    % Reconstructed: left
    plot_bp_simple(subplot(223),interp_avg_left(iF).azq_avg, ...
                   interp_avg_left(iF).elq_avg,interp_avg_left(iF).vq_norm_avg,'eckert4');
    title('Reconstructed bp: left');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    % Reconstructed: right
    plot_bp_simple(subplot(224),interp_avg_right(iF).azq_avg, ...
                   interp_avg_right(iF).elq_avg,interp_avg_right(iF).vq_norm_avg,'eckert4');
    title('Reconstructed bp: right');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    % Save figure and mat file
    figure(fig_bp_cmp)
    sfname = sprintf('%s_avg_bp_%dkHz',save_fname,D.freq.all(iF)/1e3);
    saveSameSize(fig_bp_cmp,'file',fullfile(save_path,[sfname,'.png']),...
                 'format','png','renderer','painters');
    epswrite(fullfile(save_path,[sfname,'.eps']));

end