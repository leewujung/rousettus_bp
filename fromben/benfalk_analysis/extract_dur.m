function[onsets,offsets,voc_t,I]=extract_dur(waveform,Fs,voc_t,trial_start,trial_end,pretrig_t,b2mD,manual,DIAG)
I=[];

%remove extraneous sounds below 12k
[b,a] = butter(6,30e3/(Fs/2),'high');
ddf=filtfilt(b,a,waveform);
% freqz(b,a,SR/2,SR);
data_square = smooth((ddf.^2),100);

%for marking the end time 
%we assume it's below 30k, removes some energy from echoes
[low_b, low_a]=butter(6,80e3/(Fs/2),'low'); 
waveform_low=filtfilt(low_b,low_a,ddf);
data_square_low=smooth((waveform_low.^2),100);

%for marking the start time
%we assume it's above 30k, removes some energy from previous vocs
[low_b, low_a]=butter(6,20e3/(Fs/2),'high'); 
waveform_high=filtfilt(low_b,low_a,waveform);
data_square_high=smooth((waveform_high.^2),100);

noise_length = .001*Fs; %length of data for estimating noise (1ms)
% noise_low=median(max(reshape(data_square_low(1:floor(length(data_square_low)...
%   /noise_length)*noise_length),noise_length,[])));
% noise_high=median(max(reshape(data_square_high(1:floor(length(data_square_high)...
%   /noise_length)*noise_length),noise_length,[])));

data_square_diff_high = abs(smooth(diff(data_square_high),50));
noise_diff_high=...
  median(max(reshape(data_square_diff_high(1:floor(length(data_square_diff_high)...
  /noise_length)*noise_length),noise_length,[])));
data_square_diff_low = abs(smooth(diff(data_square_low),50));
noise_diff_low=...
  median(max(reshape(data_square_diff_low(1:floor(length(data_square_diff_low)...
  /noise_length)*noise_length),noise_length,[])));

used_vocs = voc_t > trial_start & voc_t < trial_end & [1; diff(voc_t)>10e-3];
voc_t=voc_t(used_vocs);
%it's possible that you don't have a bat position for each vocalization
%even after this if you have more than one d3 trial or some ignore segments
if ~isempty(b2mD)
  b2mD=b2mD(used_vocs);
  lvl_diff = 20*log10(b2mD/.1);
end

voc_samps = round((voc_t+pretrig_t).*Fs);
buff_past = .0055*Fs;
buff_forw = .007*Fs;

if ~manual && ...
    exist('duration_threshold_fit.mat','file')
  load('duration_threshold_fit.mat');
end

onsets=nan(length(voc_t),1);
offsets=nan(length(voc_t),1);
I=nan(length(voc_t),1);
pk=nan(length(voc_t),1);
thresh1=nan(length(voc_t),1);
thresh2=nan(length(voc_t),1);
for j=1:length(voc_samps)
  voc_samp = voc_samps(j);
  first_frame=max(1,voc_samp - buff_past);
  last_frame=min(voc_samp + buff_forw,length(ddf));
  voc = ddf(first_frame:last_frame);
  data_square_voc = data_square(first_frame:last_frame);
  data_square_high_voc = data_square_high(first_frame:last_frame);
  data_square_low_voc = data_square_low(first_frame:last_frame);
  smooth_der_voc_high=data_square_diff_high(first_frame + 1:min(last_frame,length(data_square_diff_high)));
  smooth_der_voc_low=data_square_diff_low(first_frame + 1:min(last_frame,length(data_square_diff_high)));
  
  %1 ms before and after voc_t, look for maximum...
  [pk(j),loc]=max(data_square_voc( max(1,buff_past - .001*Fs):...
    min(length(data_square_voc),buff_past-.001*Fs) + .001*Fs));
  loc = loc + buff_past - .001*Fs;
  
  %no idea if this equation is correct...
  if ~isempty(b2mD)
    I(j) = 20*log10(sqrt(pk(j))) + lvl_diff(j);
  end
  
  voc_s=nan; voc_e=nan;
  
  thresh_mean_indx = -4:-1;
  thresh_mean_indx=thresh_mean_indx(j-4:j-1 >= 1);
  if j > 1
    thresh1(j)=nanmean([.9*polyval(pk_vs_thresh1.coeff,pk(j));...
      thresh1(j+thresh_mean_indx)]);
    thresh2(j)=nanmean([.9*polyval(pk_vs_thresh2.coeff,pk(j));...
      thresh1(j+thresh_mean_indx)]);
  else
    thresh1(j)=.9*polyval(pk_vs_thresh1.coeff,pk(j));
    thresh2(j)=.9*polyval(pk_vs_thresh2.coeff,pk(j));
  end
  thresh_high_noise = noise_diff_high*10;
  thresh_low_noise = noise_diff_low*10;
  
  if DIAG || manual
    figure(1); clf;
    
    hh(1)=subplot(3,1,1); cla;
    plot((1:length(voc))./Fs,voc);
    axis tight; 
    aa=axis;
    
    hh(2)=subplot(3,1,2); cla;
    [~,F,T,P] = spectrogram(voc,128,120,512,Fs);
    imagesc(T,F,10*log10(P)); set(gca,'YDir','normal');
    set(gca,'clim',[-95 -45]);
    axis tight;
    aaa=axis;
    axis([aa(1:2) aaa(3:4)]);
    
    Mhigh=max(data_square_high_voc(1:loc));
    Mlow=max(data_square_low_voc(loc:end));
    hh(3)=subplot(3,1,3); cla;
    plot((1:loc)./Fs, data_square_high_voc(1:loc)./Mhigh);
    hold on;
    plot((loc:length(data_square_low_voc))./Fs,data_square_low_voc(loc:end)...
      ./Mlow);
    th1_a=plot([0,length(data_square_high_voc)/2]./Fs,...
      [thresh1(j)./Mhigh thresh1(j)./Mhigh],'g');
    th2_a=plot([length(data_square_low_voc)/2,length(data_square_low_voc)]./Fs,...
      [thresh2(j)./Mlow thresh2(j)./Mlow],'g');
    axis tight;
    
    linkaxes(hh,'x');
    
    Mdiffhigh=max(smooth_der_voc_high(1:loc));
    Mdifflow=max(smooth_der_voc_low(loc:end));
    figure(2); clf;
    plot((1:loc)./Fs,...
      smooth_der_voc_high(1:loc)./Mdiffhigh);
    hold on;
    plot((loc:length(smooth_der_voc_low))./Fs,...
      smooth_der_voc_low(loc:end)./Mdifflow);
    axis tight;
    th12_aa=plot([1 loc loc length(smooth_der_voc_low)]./Fs,...
      [thresh_high_noise thresh_high_noise...
      thresh_low_noise thresh_low_noise]./...
      [Mdiffhigh Mdiffhigh Mdifflow Mdifflow],'g');
  end
  
  if ~manual && pk(j) > .001
    try
    [voc_s,thresh1(j), thresh_high_noise]=...
      find_thresh_crossing(data_square_high_voc(1:loc),...
      thresh1(j),0,'last',smooth_der_voc_high,thresh_high_noise);
    [voc_e,thresh2(j), thresh_low_noise]=...
      find_thresh_crossing(data_square_low_voc(loc:end),...
      thresh2(j),loc,'first',smooth_der_voc_low,thresh_low_noise);
    catch
      disp(num2str(j))
    end
  end
  
  if DIAG || manual
    figure(1);
    axes(hh(1));
    hold on;
    plot([voc_s voc_s]./Fs,[min(voc) max(voc)],'r');
    plot([voc_e voc_e]./Fs,[min(voc) max(voc)],'r');
    axis tight;
    aa=axis;
    
    axes(hh(2));
    hold on;
    plot([voc_s voc_s]./Fs,[0 Fs/2],'r');
    plot([voc_e voc_e]./Fs,[0 Fs/2],'r');
    axis tight;
    aaa=axis;
    axis([aa(1:2) aaa(3:4)]);
    
    axes(hh(3));
    plot([voc_s voc_s]./Fs,...
      [0 1],'r');
    plot([voc_e voc_e]./Fs,...
      [0 1],'r');
    delete([th1_a th2_a]);
    plot([0,length(data_square_high_voc)/2]./Fs,...
      [thresh1(j)./Mhigh thresh1(j)./Mhigh],'g');
    plot([length(data_square_low_voc)/2,length(data_square_low_voc)]./Fs,...
      [thresh2(j)./Mlow thresh2(j)./Mlow],'g');
    axis tight;
    
    Mdiffhigh=max(smooth_der_voc_high(1:loc));
    Mdifflow=max(smooth_der_voc_low(loc:end));
    figure(2);
    axis tight;
    aa = axis;
    plot([voc_s voc_s]./Fs,[0 aa(4)],'r');
    plot([voc_e voc_e]./Fs,[0 aa(4)],'r');
    delete(th12_aa);
    plot([1 loc loc length(smooth_der_voc_low)]./Fs,...
      [thresh_high_noise thresh_high_noise...
      thresh_low_noise thresh_low_noise]./...
      [Mdiffhigh Mdiffhigh Mdifflow Mdifflow],'g');
  
%     figure(3), clf;
%     plot((1:length(data_square_low_voc))./Fs,data_square_low_voc);
%     hold on;
% %     plot([voc_s voc_s]./Fs,[0 max(data_square_low)],'r');
%     plot([voc_e voc_e]./Fs,[0 max(data_square_low_voc)],'r');
%     plot([0,length(data_square_low_voc)]./Fs,[thresh2 thresh2],'g');
%     axis tight;
    
%     figure(4), hold on;
%     scatter(pk,polyval(thresh1fit.coeff,pk));
    
%     figure(3); clf;
%     plot(data_square_low)
    
    %     if ~manual && ~isnan(voc_s) && ~isnan(voc_e)
    %       figure(4); clf;
    %       NFFT=1024;
    %       Y=fft(voc(voc_s:voc_e),NFFT);
    %       f=Fs/2*linspace(0,1,NFFT/2+1);
    %       plot(f,smooth(2*abs(Y(1:NFFT/2+1)),10))
    %       axis tight;
    %     end
    
    %     figure(5); clf;
    %     scatter(b2mD,sqrt(data_square(voc_samps)));
    %     hold on;
    %     plot(b2mD(k),sqrt(data_square(voc_samps(k))),'or','markerfacecolor','r');
    
    %     figure(6); clf;
    %     scatter(I,offsets-onsets);
    
    if manual
      disp('mark start and stop of vocalization')
      axes(hh(3));
      [x]=ginput(2);
      if diff(x)<0
        disp('wrong order')
        x=ginput(2);
      end
      voc_s = x(1)*Fs;
      thresh1(j)=data_square_high_voc(x(1));
      voc_e = x(2)*Fs;
      thresh2(j)=data_square_low_voc(x(2));
    end
  end
  
  onsets(j)=(voc_s + voc_samp - buff_past)/Fs - pretrig_t;
  offsets(j)=(voc_e + voc_samp - buff_past)/Fs - pretrig_t;
end


function[cross,thresh,noise_thresh]=find_thresh_crossing(data,thresh,offset,type,data_diff,noise_thresh)
if strcmp(type,'last')
  init_offset = max(1,length(data)-500);
  ii=find(data(init_offset:end)>max(data(init_offset:end))*.75,1,type)+init_offset-1;
  M=data(ii);
  thresh = mean([abs(thresh) median(data)*1.5]);
  cross=find(data(1:ii) < thresh,1,type);
  while isempty(cross)
    thresh = thresh*1.5;
    cross = find(data(1:ii) < thresh,1,type);
    if thresh > M
      break;
    end
  end
  [~,locs]=findpeaks(data(1:cross),'MINPEAKHEIGHT',max(data(1:cross))/2,...
    'NPEAKS',3,'THRESHOLD',max(data(1:cross))/1e3);
  if ~isempty(locs) && (cross - locs(end)) < 250
    [cross, thresh]=find_thresh_crossing(data(1:locs(end)),thresh,0,type,...
      data_diff,noise_thresh);
  end
  cross_diff=find(data_diff(1:cross) < noise_thresh,1,type);
  if ~isempty(cross_diff)
    cross=cross_diff;
  end
else
  ii=find(data>max(data(1:min(500,length(data))))*.75,1,type);
  M=data(ii);
  thresh = min([thresh median(data)*2]);
  cross=find(data(ii:end) < thresh,1,type)+ii-1;
  while isempty(cross)
    thresh = thresh*1.5;
    cross = find(data(ii:end) < thresh,1,type)+ii-1;
    if thresh > M
      break;
    end
  end
  noise_thresh=min([noise_thresh median(data_diff)*2]);
  strt_indx = max([0 cross-85 ii])+offset;
  cross=find(data_diff(strt_indx:end) < noise_thresh,1,type)...
    + strt_indx - offset;
end

cross = cross + offset;




