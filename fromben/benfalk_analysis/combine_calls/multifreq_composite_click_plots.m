% 2015 12 28  Plots for multi-freq composite clicks

clear

usrn = getenv('username');
addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);

save_plot_opt = 1;

% Load compiled rotated data
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
composite_path = '20151228_composite_click_figs_3bats';
compile_file = 'rousettus_36134_composite_data.mat';
load(fullfile(base_path,composite_path,compile_file));

% Path/filename for saving data and figures
save_header = 'rousettus_36134';
save_path = fullfile(base_path,'20160209_composite_multifreq');  % set path for saving files
if ~exist(save_path,'dir')
    mkdir(save_path);
end

num_freq = length(param.freq_wanted);

%% Scatter plot for all frequencies to show raw data
for iF=1:num_freq
    fig_scat = figure;
    set(fig_scat,'position',[70,100,1200,560]);
    
    subplot(121)  % left click
    scatter(scatter_data.left.az_shift_tilt,...
            scatter_data.left.el_shift_tilt,15,...
            scatter_data.left.call_dB_norm(iF,:),'fill')
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    title(sprintf('All left click data points, %d kHz',param.freq_wanted(iF)/1e3));
    colorbar('location','southoutside');
    grid
    
    subplot(122)  % right click
    scatter(scatter_data.right.az_shift_tilt,...
            scatter_data.right.el_shift_tilt,15,...
            scatter_data.right.call_dB_norm(iF,:),'fill')
    xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
    title(sprintf('All right click data points, %d kHz',param.freq_wanted(iF)/1e3));
    colorbar('location','southoutside');
    grid
    
    if save_plot_opt == 1
        saveSameSize(fig_scat,'file',fullfile(save_path,[save_header,'_',sprintf('all_pts_freq%dkHz.png',param.freq_wanted(iF)/1e3)]),...
            'format','png','renderer','painters');
    end
end



%% Averaged composite beampattern for all frequencies
for iF=1:num_freq
    fig_comp_avg = figure;
    set(fig_comp_avg,'position',[60,190,1250,390]);
    
    subplot(121)  % left click
    axesm eckert4
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
    contourfm(averaged_composite.left.interp(iF).elq_avg,...
              averaged_composite.left.interp(iF).azq_avg,...
              averaged_composite.left.interp(iF).vq_norm_avg,-30:3:-3);
    caxis([-30 0])
    colormap(parula(10))
    colorbar('Ticks',-30:3:0,'location','southoutside');
    title(sprintf('Averaged left clicks, %d kHz, th=%d, bin=%ddeg',...
                  param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));

    subplot(122)  % right click
    axesm eckert4
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
    contourfm(averaged_composite.right.interp(iF).elq_avg,...
              averaged_composite.right.interp(iF).azq_avg,...
              averaged_composite.right.interp(iF).vq_norm_avg,-3:-3:-30);
    caxis([-30 0])
    colormap(parula(10))
    colorbar('Ticks',-30:3:0,'location','southoutside');
    title(sprintf('Averaged right clicks, %d kHz, th=%d, bin=%ddeg',...
                   param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));
               
    if save_plot_opt == 1
        saveSameSize(fig_comp_avg,'file',fullfile(save_path,[save_header,'_',sprintf('composite_avg_freq%dkHz.png',param.freq_wanted(iF)/1e3)]),...
            'format','png','renderer','painters');
    end
end


%% Raw composite beampattern for all frequencies
for iF=1:num_freq
    fig_comp_raw = figure;
    set(fig_comp_raw,'position',[60,190,1250,390]);
    
    subplot(121)  % left click
    axesm eckert4
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
%     contourfm(raw_composite.left.el/pi*180,...
%               raw_composite.left.az/pi*180,...
%               squeeze(raw_composite.left.vq(iF,:,:)),-30:3:-3);
    contourfm(raw_composite.left.elq/pi*180,...
              raw_composite.left.azq/pi*180,...
              squeeze(raw_composite.left.vq(iF,:,:)),-30:3:-3);
    caxis([-30 0])
    colormap(parula(10))
    colorbar('Ticks',-30:3:0,'location','southoutside');
    title(sprintf('Merged left clicks, %d kHz, no averaging',param.freq_wanted(iF)/1e3));

    subplot(122)  % right click
    axesm eckert4
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    axis off
%     contourfm(raw_composite.right.el/pi*180,...
%               raw_composite.right.az/pi*180,...
%               squeeze(raw_composite.right.vq(iF,:,:)),-3:-3:-30);
    contourfm(raw_composite.right.elq/pi*180,...
              raw_composite.right.azq/pi*180,...
              squeeze(raw_composite.right.vq(iF,:,:)),-3:-3:-30);
    caxis([-30 0])
    colormap(parula(10))
    colorbar('Ticks',-30:3:0,'location','southoutside');
    title(sprintf('Merged right clicks, %d kHz, no averaging',param.freq_wanted(iF)/1e3));
    
    if save_plot_opt == 1
        saveSameSize(fig_comp_raw,'file',fullfile(save_path,[save_header,'_',sprintf('composite_raw_freq%dkHz.png',param.freq_wanted(iF)/1e3)]),...
            'format','png','renderer','painters');
    end
end


%% Plot multi-freq contour
num_freq = length(param.freq_wanted);
colorset = jet(num_freq);

fig_multifreq_contour = figure;
set(fig_multifreq_contour,'position',[150 150 1200 400]);
for iF=1:num_freq
    subplot(121)  % left click
    hold on
    plot(multifreq_3dB_contour.left{iF}(:,1),multifreq_3dB_contour.left{iF}(:,2),'linewidth',2,'color',colorset(iF,:));
    subplot(122)  % right click
    hold on
    plot(multifreq_3dB_contour.right{iF}(:,1),multifreq_3dB_contour.right{iF}(:,2),'linewidth',2,'color',colorset(iF,:));
end
subplot(121)
axis equal
colormap(jet(num_freq))
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
         'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
axis([param.map.map_plot_xlim,param.map.map_plot_ylim])
grid
subplot(122)
axis equal
colormap(jet(num_freq))
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
         'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
axis([param.map.map_plot_xlim,param.map.map_plot_ylim])
grid

if save_plot_opt==1
    saveSameSize(fig_multifreq_contour,'file',fullfile(save_path,[save_header,'_','multifreq_3dB.png']),...
        'format','png','renderer','painters');
end

