clear; close all;

call_type='eptesicus';
% call_type='rousettus';

if strcmp(call_type,'eptesicus')
  bp_path='..\proc_output_eptesicus_new\';
  load([bp_path 'eptesicus_20150824_LB62_09_mic_data_bp_proc.mat'])

  mic_data_path='..\mic_data\';
  load([mic_data_path 'eptesicus_20150824_LB62_09_mic_data.mat'])
  callnum=15;
else
  bp_path='..\proc_output\';
  load([bp_path 'rousettus_20150825_36134_02_mic_data_bp_proc.mat'])

  mic_data_path='..\mic_data\';
  load([mic_data_path 'rousettus_20150825_36134_02_mic_data.mat'])
  % callnum=19;
  callnum=19;
end

t=((1:length(sig))-1)/fs-4;
callt=proc.call_receive_time(callnum);

% if strcmp(call_type,'rousettus')
%   idx=t>callt-.0021 & t<callt+.0012;
% else
  idx=t>callt-.001 & t<callt+.0023;
% end

% [~,best_ch]=max(max(sig(idx,:)));
best_ch=33; %the 1/8 inch calib. microphone
voc=sig(idx,best_ch);

colordef white

figure(1); clf; set(gcf,'pos',[10 40 400 700],'color','white')
ax1=subaxis(3,1,1,'m',.02,'sv',.01,'ml',.2,'mb',.1);
plot(((1:length(voc))-1)/fs*1e3,voc);
YL(1)=ylabel('Intensity');
set(gca,'fontsize',14)

ax2=subaxis(3,1,2,'m',.02,'sv',.01,'ml',.2,'mb',.1);
[~,F,T,P] = spectrogram(voc,64,60,1e3,fs);
imagesc(T*1e3,F,10*log10(P)); 
set(gca,'YDir','normal');
ca=caxis;
caxis([ca(2)-25,ca(2)])
linkaxes([ax1 ax2],'x')
axis tight
% a=axis();
% axis([ a(1:3) 100e3])
colormap hot
cmap = colormap;
colormap(flipud(cmap))
xlabel('Time (ms)')
YL(2)=ylabel('Frequency (kHz)');
set(gca,'YTick',(25:25:125) *1e3,'yticklabel',(25:25:125),'fontsize',14)

subaxis(3,1,3,'m',.02,'sv',.21,'ml',.2,'mb',.1);
% pwelch(voc,[],[],[],fs)
% axis tight;
% a=axis;
% axis([a(1) 100 a(4)-30 a(4)])
% ylabel('Power/freq (dB/Hz)')
% YL(3)=get(gca,'YLabel');
% title('Power Spectral Density','fontweight','normal')

cla
spec=proc.call_psd_dB_comp_re20uPa_withbp{callnum,best_ch};
freq=proc.call_freq_vec{callnum,best_ch};
plot(freq./1e3,spec);
axis tight;
% a=axis;
% axis([a(1) 100e3 a(4)-30 a(4)])
ylabel('dB')
xlabel('Freq (kHz)')
YL(3)=get(gca,'YLabel');
title('Spectrum','fontweight','normal')
set(gca,'fontsize',14)


for k=1:3
  YL(k).Units = 'normalized';
  ypos(k,:)=get(YL(k),'position');
end
for k=1:3
  YL(k).Position = [min(ypos(:,1)) ypos(k,2:3)];
end


