% 2016 08 08  Assemble composite clicks from simulated mic receptions
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 09 23  Enable assembling composite click from multi-freq simulation

clear

usrn = getenv('username');
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
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
A.map = D.map;

% Set save_fname
save_fname = sprintf('%s_%s_bin%d_th%d',...
    script_name,noise_std,binsize,threshold);

% Merge all data
click_side_all = nan(length(simu_file),num_ch);
az_elpctr_tilt_all = nan(length(simu_file),num_ch);
el_elpctr_tilt_all = nan(length(simu_file),num_ch);
call_dB_all = nan(length(simu_file),num_ch,length(D.freq.all));
for iS=1:50%length(simu_file)

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
call_dB_merge = reshape(call_dB_all,[],num_freq);

idx_right_good = find(click_side_merge==1 & ~isnan(call_dB_merge(:,1)));
call_dB_right = call_dB_merge(idx_right_good,:);
az_right = az_elpctr_tilt_all(idx_right_good);
el_right = el_elpctr_tilt_all(idx_right_good);

idx_left_good = find(click_side_merge==0 & ~isnan(call_dB_merge(:,1)));
call_dB_left = call_dB_merge(idx_left_good,:);
az_left = az_elpctr_tilt_all(idx_left_good);
el_left = el_elpctr_tilt_all(idx_left_good);



% Composite multi-freq model click avg_bp and -3dB contours
for iF=1:num_freq
    fprintf('Merging data at freq %d kHz\n',D.freq.all(iF)/1e3);

    % Interpolation using averaged data
    [interp_avg_right(iF),bin_avg_right(iF)] = ...
        average_call(az_right,el_right,call_dB_right(:,iF),binsize,'eckert4',threshold);
    [interp_avg_left(iF),bin_avg_left(iF)] = ...
        average_call(az_left,el_left,call_dB_left(:,iF),binsize,'eckert4',threshold);

    % Find multi-freq -3dB contours for composite model clicks
    [azq,elq,vq] = get_ortho_grid_azel(interp_avg_left(iF).azq_avg,...
                                       interp_avg_left(iF).elq_avg,...
                                       interp_avg_left(iF).vq_norm_avg);
    [~,c_level_nan_left] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
    [c3db_xy_left(:,1),c3db_xy_left(:,2)] = ...
        mfwdtran(D.map.mstruct,c_level_nan_left(:,2),c_level_nan_left(:,1));  % [az,el] to [x,y]

    [azq,elq,vq] = get_ortho_grid_azel(interp_avg_right(iF).azq_avg,...
                                       interp_avg_right(iF).elq_avg,...
                                       interp_avg_right(iF).vq_norm_avg);
    [~,c_level_nan_right] = get_main_contour(vq,unique(azq(:)),unique(elq(:)),-3);
    [c3db_xy_right(:,1),c3db_xy_right(:,2)] = ...
        mfwdtran(D.map.mstruct,c_level_nan_right(:,2),c_level_nan_right(:,1));  % [az,el] to [x,y]

    c3db_xy_all_right{iF} = c3db_xy_right;
    c3db_xy_all_left{iF} = c3db_xy_left;

    clear c3db_xy_left c3db_xy_right
end

% Save to struct
A.averaged_composite.right.interp = interp_avg_right;
A.averaged_composite.right.bin = bin_avg_right;
A.averaged_composite.left.interp = interp_avg_left;
A.averaged_composite.left.bin = bin_avg_left;

A.multifreq_3dB_contour.right = c3db_xy_all_right;
A.multifreq_3dB_contour.left = c3db_xy_all_left;


% Find beam center for composite model clicks
for iF=1:num_freq
    fprintf('Finding beam center at freq %d kHz\n',D.freq.all(iF)/1e3);
    
    % Find beam center based on energy
    % --- left composite click
    xx = interp_avg_left(iF).vq_norm_avg(:);
    [left_max_val,left_max_idx] = max(xx);
    left(iF).max.el = interp_avg_left(iF).elq_avg(left_max_idx);
    left(iF).max.az = interp_avg_left(iF).azq_avg(left_max_idx);
    xx(isnan(xx)) = -inf;
    [~,left_sort_idx] = sort(xx,'descend');
    ii = xx(left_sort_idx)>-1;
    left(iF).top.el = mean(interp_avg_left(iF).elq_avg(left_sort_idx(ii)));
    left(iF).top.az = mean(interp_avg_left(iF).azq_avg(left_sort_idx(ii)));
    % --- right composite click
    xx = interp_avg_right(iF).vq_norm_avg(:);
    [right_max_val,right_max_idx] = max(xx);
    right(iF).max.el = interp_avg_right(iF).elq_avg(right_max_idx);
    right(iF).max.az = interp_avg_right(iF).azq_avg(right_max_idx);
    xx(isnan(xx)) = -inf;
    [~,right_sort_idx] = sort(xx,'descend');
    ii = xx(right_sort_idx)>-1;
    right(iF).top.el = mean(interp_avg_right(iF).elq_avg(right_sort_idx(ii)));
    right(iF).top.az = mean(interp_avg_right(iF).azq_avg(right_sort_idx(ii)));

    % Fit ellipse
    % --- left composite click
    [left_raw,left_rot_max,left_rot_elpctr,left_rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(bin_avg_left(iF).az_avg,...
                                  bin_avg_left(iF).el_avg,...
                                  bin_avg_left(iF).avg_call_dB,...
                                  D.map.map_projection,0.005);
    [left_el_ectr,left_az_ectr] = minvtran(D.map.mstruct,left_rot_max.E.x0, ...
                                           left_rot_max.E.y0);  % inverse map projection

    [left(iF).ectr.el,left(iF).ectr.az] = rotatem(left_el_ectr,left_az_ectr,...
                                    [left(iF).max.el,left(iF).max.az],...
                                    'inverse','degrees');
    % --- right composite click
    [right_raw,right_rot_max,right_rot_elpctr,right_rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(bin_avg_right(iF).az_avg,...
                                  bin_avg_right(iF).el_avg,...
                                  bin_avg_right(iF).avg_call_dB,...
                                  D.map.map_projection,0.005);
    [right_el_ectr,right_az_ectr] = minvtran(D.map.mstruct,right_rot_max.E.x0,...
                                             right_rot_max.E.y0);  % inverse map projection
    [right(iF).ectr.el,right(iF).ectr.az] = rotatem(right_el_ectr,right_az_ectr,...
                                    [right(iF).max.el,right(iF).max.az],...
                                    'inverse','degrees');
end
A.bpctr.right = right;
A.bpctr.left = left;



% Save mat file
save(fullfile(save_path,[save_fname,'.mat']),'-struct','A');



% Plot params
suptitle_str_pre = sprintf('%s %s',model_shape,noise_std);
cgrey = 200*ones(1,3)/255;
cvec=-30:3:0;
num_freq_plot = num_freq-2;
colorset = jet(num_freq_plot);  % colormap used for contours
contour_sm_len = 10;


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



% Plot multi-freq -3dB contours
bpctr_opt = 'ectr';

fig_left = figure;
axesm eckert4
axis off
framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
gridm('gcolor',cgrey,'glinestyle','-');
hold on

fig_right = figure;
axesm eckert4
axis off
framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
gridm('gcolor',cgrey,'glinestyle','-');
hold on

[~,I] = sort(D.freq.all);
cnt = 0;
for iF = I(2:end-1)  % from 25 to 55 kHz
    cnt = cnt+1;

    % Get contours
    xy_sm_left(:,1) = smooth(c3db_xy_all_left{iF}(:,1),contour_sm_len);
    xy_sm_left(:,2) = smooth(c3db_xy_all_left{iF}(:,2),contour_sm_len);
    xy_sm_left(isnan(xy_sm_left(:,1)),:) = NaN;
    
    xy_sm_right(:,1) = smooth(c3db_xy_all_right{iF}(:,1),contour_sm_len);
    xy_sm_right(:,2) = smooth(c3db_xy_all_right{iF}(:,2),contour_sm_len);
    xy_sm_right(isnan(xy_sm_right(:,1)),:) = NaN;

    % Left click
    figure(fig_left)
    % --- -3dB contour
    plot(xy_sm_left(:,1),xy_sm_left(:,2),'linewidth',3,'color',colorset(cnt,:));
    switch bpctr_opt
      case 'max'  % Plot max beam energy point
        plotm(left(iF).max.el,left(iF).max.az,'x','markersize', ...
              8,'linewidth',2,'color',colorset(cnt,:));
      case 'top'  % Plot dB>-1 beam energy point
        plotm(left(iF).top.el,left(iF).top.az,'^','markersize', ...
              8,'linewidth',2,'color',colorset(cnt,:));
      case 'ectr' % location of center of best-fitting ellipse
        plotm(left(iF).ectr.el,left(iF).ectr.az,'o', ...
              'markersize',8,'linewidth',2,'color',colorset(cnt,:));
    end

    % Right click
    figure(fig_right)
    % --- -3dB contour
    plot(xy_sm_right(:,1),xy_sm_right(:,2),'linewidth',3,'color',colorset(cnt,:));
    switch bpctr_opt
      case 'max'  % Plot max beam energy point
        plotm(right(iF).max.el,right(iF).max.az,'x','markersize', ...
              8,'linewidth',2,'color',colorset(cnt,:));
      case 'top'  % Plot dB>-1 beam energy point
        plotm(right(iF).top.el,right(iF).top.az,'^','markersize', ...
              8,'linewidth',2,'color',colorset(cnt,:));
      case 'ectr' % location of center of best-fitting ellipse
        plotm(right(iF).ectr.el,right(iF).ectr.az,'o', ...
              'markersize',8,'linewidth',2,'color',colorset(cnt,:));
    end
    
    clear xy_sm_*
end


figure(fig_left)
colormap(jet(num_freq_plot))
colorbar('Ticks',linspace(0+1/num_freq_plot/2,1-1/num_freq_plot/2,num_freq_plot),...
    'TickLabels',{num2str(D.freq.all(I(2:end-1))'/1e3)},'location','southoutside');
grid
tightmap
title('Averaged left click');

sfname = sprintf('%s_cntr_left_%s',save_fname,bpctr_opt);
saveas(fig_left,fullfile(save_path,[sfname,'.fig']),'fig');
saveSameSize_res(fig_left,150,'file',fullfile(save_path,[sfname,'.png']),...
                 'format','png','renderer','painters');
epswrite(fullfile(save_path,[sfname,'.eps']));


figure(fig_right)
colormap(jet(num_freq_plot))
colorbar('Ticks',linspace(0+1/num_freq_plot/2,1-1/num_freq_plot/2,num_freq_plot),...
    'TickLabels',{num2str(D.freq.all(I(2:end-1))'/1e3)},'location','southoutside');
grid
tightmap
title('Averaged right click');

sfname = sprintf('%s_cntr_right_%s',save_fname,bpctr_opt);
saveas(fig_right,fullfile(save_path,[sfname,'.fig']),'fig');
saveSameSize_res(fig_right,150,'file',fullfile(save_path,[sfname,'.png']),...
                 'format','png','renderer','painters');
epswrite(fullfile(save_path,[sfname,'.eps']));
