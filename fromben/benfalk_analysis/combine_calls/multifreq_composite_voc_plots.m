% 2015 12 28  Plots for multi-freq composite clicks

clear, close all;

% bat='39184';
% bat='34271';
bat='36134';
species='rousettus';

% bat='LB88';
% bat='LB62';
% bat='LB53';
% species='eptesicus';

addpath('E:\repositories\beampattern_processing\fit_beams\');
addpath('..\')

save_plot_opt = 1;
scatter_all_freq=0;
avg_composite_bp=0;
raw_composite_bp=0;
multi_freq_contour=1;


% Load compiled rotated data
base_path = '..\..\';
composite_path = 'composite_figs';
compile_file = [bat '_composite_data.mat'];
load(fullfile(base_path,composite_path,compile_file));

% Path/filename for saving data and figures
save_header = [species '_' bat];
save_path = fullfile(base_path,'composite_figs');  % set path for saving files
if ~exist(save_path,'dir')
  mkdir(save_path);
end

num_freq = length(param.freq_wanted);

%% Scatter plot for all frequencies to show raw data
if scatter_all_freq
  for iF=1:num_freq
    fig_scat = figure;
    set(fig_scat,'position',[10,50,650,650]);
    
    if strcmp(species,'rousettus')
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
    else
      scatter(scatter_data.all.az_shift_tilt,...
        scatter_data.all.el_shift_tilt,15,...
        scatter_data.all.call_dB_norm(iF,:),'fill')
      caxis([-30,0])
      xlabel('Azimuth (deg)'); ylabel('Elevation (deg)');
      title(sprintf('All data points, %d kHz',param.freq_wanted(iF)/1e3));
      colorbar('location','southoutside');
      grid on, axis equal;
      axis([-90 90 -90 90]);
    end
    drawnow;
    if save_plot_opt
      export_fig(fullfile(save_path,[save_header,'_',...
        sprintf('all_pts_%dkHz.png',param.freq_wanted(iF)/1e3)]),...
        '-png','-nocrop');
      %     saveSameSize(fig_scat,'file',fullfile(save_path,[save_header,'_',sprintf('all_pts_freq%dkHz.png',param.freq_wanted(iF)/1e3)]),...
      %       'format','png','renderer','painters');
    end
  end
  close all;
end

%% Averaged composite beampattern for all frequencies
if avg_composite_bp
  for iF=1:num_freq
    fig_comp_avg = figure;
    set(fig_comp_avg,'position',[10,50,650,650]);
    
    if strcmp(species,'rousettus')
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
    else
      cla
      axesm eckert4
      framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
      gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
      axis off
      contourfm(averaged_composite.all.interp(iF).elq_avg,...
        averaged_composite.all.interp(iF).azq_avg,...
        averaged_composite.all.interp(iF).vq_norm_avg,-3:-3:-30);
      colormap(parula(10))
      caxis([-30 0])
      colorbar('Ticks',-30:3:0,'location','southoutside');
      axis([-pi/2 pi/2 -pi/2 pi/2])
      title(sprintf('Averaged, %d kHz, th=%d, bin=%ddeg',...
        param.freq_wanted(iF)/1e3,param.composite_threshold,param.composite_binsize));
    end
    
    drawnow
    if save_plot_opt
      export_fig(fullfile(save_path,[save_header,'_',...
        sprintf('composite_avg_%dkHz.png',param.freq_wanted(iF)/1e3)]),...
        '-png','-nocrop');
      %     saveSameSize(fig_comp_avg,'file',fullfile(save_path,[save_header,'_',...
      %       sprintf('composite_avg_freq%dkHz.png',param.freq_wanted(iF)/1e3)]),...
      %       'format','png','renderer','painters');
    end
  end
  close all;
end

%% Raw composite beampattern for all frequencies
if raw_composite_bp
  for iF=1:num_freq
    fig_comp_raw = figure;
    set(fig_comp_raw,'position',[10,50,650,650]);
    
    if strcmp(species,'rousettus')
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
      title(sprintf('Merged left clicks, %d kHz, no averaging',...
        param.freq_wanted(iF)/1e3));
      
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
      title(sprintf('Merged right clicks, %d kHz, no averaging',...
        param.freq_wanted(iF)/1e3));
    else
      axesm eckert4
      framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
      gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
      axis off
      %     contourfm(raw_composite.all.el/pi*180,...
      %               raw_composite.all.az/pi*180,...
      %               squeeze(raw_composite.all.vq(iF,:,:)),-3:-3:-30);
      contourfm(raw_composite.all.elq/pi*180,...
        raw_composite.all.azq/pi*180,...
        squeeze(raw_composite.all.vq(iF,:,:)),-3:-3:-30);
      caxis([-30 0])
      colormap(parula(10))
      colorbar('Ticks',-30:3:0,'location','southoutside');
      axis([-pi/2 pi/2 -pi/2 pi/2])
      title(sprintf('Merged, %d kHz, no averaging',...
        param.freq_wanted(iF)/1e3));
    end
    
    drawnow
    if save_plot_opt
      export_fig(fullfile(save_path,[save_header,'_',...
        sprintf('composite_raw_%dkHz.png',param.freq_wanted(iF)/1e3)]),...
        '-png','-nocrop');
      %     saveSameSize(fig_comp_raw,'file',fullfile(save_path,[save_header,'_',...
      %       sprintf('composite_raw_freq%dkHz.png',param.freq_wanted(iF)/1e3)]),...
      %       'format','png','renderer','painters');
    end
  end
  close all;
end

%% Plot multi-freq contour
if multi_freq_contour
  num_freq = length(param.freq_wanted);
  %divergent colormap: http://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=9
  colorset = [215 48 39; 244 109 67; 253 174 97; 254 224 144; 255 255 191; ...
    224 243 248; 171 217 233; 116 173 209; 69 117 180;]./255;
  colorset = flipud(colorset);
  
  % colorset = parula(num_freq);
  
  fig_multifreq_contour = figure;
  set(fig_multifreq_contour,'position',[10 50 600 600],'color','w');
  cla; hold on;
  for iF=1:num_freq
    
    if strcmp(species,'rousettus')
      set(fig_multifreq_contour,'position',[10 50 1200 600]);
      subplot(121)  % left click
      hold on;
      plot(multifreq_3dB_contour.left{iF}(:,1)*180/pi,...
        multifreq_3dB_contour.left{iF}(:,2)*180/pi,...
        'linewidth',2,'color',colorset(iF,:));
      subplot(122)  % right click
      hold on;
      plot(multifreq_3dB_contour.right{iF}(:,1)*180/pi,...
        multifreq_3dB_contour.right{iF}(:,2)*180/pi,...
        'linewidth',2,'color',colorset(iF,:));
    else
      plot(multifreq_3dB_contour.all{iF}(:,1)*180/pi,...
        multifreq_3dB_contour.all{iF}(:,2)*180/pi,...
        'linewidth',2,'color',colorset(iF,:));
    end
    
  end
  
  if strcmp(species,'rousettus')
    subplot(121)
    axis equal
    colormap(colorset)
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
      'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
    axis([-90 90 -90 90])
    grid
    subplot(122)
    axis equal
    colormap(colorset)
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
      'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
    axis([-90 90 -90 90])
    grid
  else
    axis equal
    colormap(colorset)
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
      'TickLabels',{num2str(param.freq_wanted'/1e3)},'location','southoutside');
    axis([-90 90 -90 90])
    grid
  end
  if save_plot_opt
    export_fig(fullfile(save_path,[save_header,'_','multifreq_3dB.png']),...
      '-png','-nocrop');
    %   saveSameSize(fig_multifreq_contour,'file',fullfile(save_path,[save_header,'_',...
    %     'multifreq_3dB.png']),'format','png','renderer','painters');
  end
end
