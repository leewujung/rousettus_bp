clear; close all;

% bp_path='..\proc_output_eptesicus_new\';
% load([bp_path 'eptesicus_20150824_LB62_09_mic_data_bp_proc.mat'])
% 
% mic_data_path='..\mic_data\';
% load([mic_data_path 'eptesicus_20150824_LB62_09_mic_data.mat'])
% call_indx=15;

bp_path='..\proc_output\';
load([bp_path 'rousettus_20150825_36134_02_mic_data_bp_proc.mat'])

% mic_data_path='..\mic_data\';
% load([mic_data_path 'rousettus_20150825_36134_02_mic_data.mat'])
call_indx=19;

%divergent colormap: http://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=9
colorset = [215,48,39;252,141,89;254,224,144;255,255,191;224,243,248;145,191,219;69,117,180;]./255;
colorset = flipud(colorset);

figure(1); clf, set(gcf,'pos',[10 40 460 520])
hold on; grid on;
freqs=30:5:60;

for frq=1:length(freqs)
  freq_desired=freqs(frq);
  interp_method='rb_rbf';

  goodch=1:length(mic_loc(:,1));
  mic_num = 1:mic_data.num_ch_in_file;

  call_dB = nan(1,mic_data.num_ch_in_file);
  for iM=mic_num
    freq = proc.call_freq_vec{call_indx,iM};
    [~,fidx] = min(abs(freq-freq_desired*1e3));
    call_dB(iM) = proc.call_psd_dB_comp_re20uPa_withbp{call_indx,iM}(fidx);
  end

  ch_ex_sig = find(isnan(call_dB)); % low quality channel from call extraction function
  ch_good_loc = ~isnan(mic_loc(:,1))';  % Check for mics without location data

  angle_notnanidx = ~ismember(mic_num,ch_ex_sig) & ch_good_loc;

  mic_to_bat_angle = squeeze(proc.mic_to_bat_angle(call_indx,:,:));
  az = mic_to_bat_angle(angle_notnanidx,1);
  el = mic_to_bat_angle(angle_notnanidx,2);

  %       [azq,elq] = meshgrid(linspace(min(az),max(az),150),linspace(min(el),max(el),150));
  [azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
  if strcmp(interp_method,'rb_natural')  % natural neighbor interpolation
    vq = griddata(az,el,call_dB(angle_notnanidx),azq,elq,'natural');
  elseif strcmp(interp_method,'rb_rbf')  % radial basis function interpolation
    vq = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],call_dB(angle_notnanidx),'RBFFunction','multiquadrics'));
    vq = reshape(vq,size(azq));
  end
  maxref = max(call_dB(angle_notnanidx));
  vq_norm = vq-maxref;

  % Find indices within measured polygon
  k = boundary(az,el,0);  % outer boundary of all measured points
  [in,on] = inpolygon(azq,elq,az(k),el(k));
  in_smpl_poly = in|on;
  clear in on
  vq(~in_smpl_poly) = NaN;
  vq_norm(~in_smpl_poly) = NaN;

  [~,II] = max(vq_norm(:));
  [I,J]=ind2sub(size(vq_norm),II);

  elqm = elq;
  elqm(isnan(vq_norm)) = NaN;
  azqm = azq;
  azqm(isnan(vq_norm)) = NaN;

  contour(azq/pi*180,elq/pi*180,vq_norm,'linewidth',2,'levellist',-3,...
    'linecolor',colorset(frq,:))

end

axis([-90 90 -90 90]);
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
colormap(colorset)
CC=colorbar('ticks',freqs,'location','southoutside');
ylabel(CC, 'Frequency')
caxis([freqs(1)-2.5,freqs(end)+2.5])
set(gca,'fontsize',14)




% 
% %plotting the entire beamshape
% cla
% axesm eckert4;
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
% gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% axis off
% geoshow(elqm/pi*180,azqm/pi*180,vq_norm,'displaytype','texturemap');
% contourm(elq/pi*180,azq/pi*180,vq_norm,-3,'w','linewidth',2);
% scatterm(el/pi*180,az/pi*180,[],[0 0 0],'+')
% %       textm(el/pi*180,az/pi*180,num2str(mic_num(angle_notnanidx)'),'horizontalalignment','center');
% cc = [-30 0];
% caxis(cc);
% colorbar('southoutside');
% %       tightmap
% %       setm(gca,'ParallelLabel','on','MeridianLabel','on')
% 
% vq_norm_min = min(vq_norm(:));
% contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
% cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');
