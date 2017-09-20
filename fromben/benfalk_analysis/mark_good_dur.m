function [saving,new_duration_data]=mark_good_dur(dur_data,waveform,Fs,pretrig_t,new_duration_data)
saving = 0;

disp('Type Enter/Return to keep vocalization');
disp('minus to delete voc, ESC to quit');
disp('Arrows LEFT/RIGHT to move Back/Forward');

ff=figure(1);
set(ff,'position',[30 45 650 700])

vv_indx=find(~isnan(new_duration_data(:,2)))';

max_dur=max(dur_data(vv_indx(:),3)-...
  dur_data(vv_indx(:),2));

vv=1;
buffer_s = round((2e-3).*Fs);
buffer_e = round((2e-3).*Fs);
while vv <= length(vv_indx)
  voc_s = dur_data(vv_indx(vv),2);
  voc_e = dur_data(vv_indx(vv),3);
  samp_s = round((voc_s + pretrig_t).*Fs);
  samp_e = round((voc_e + pretrig_t).*Fs);
  
  clf(ff);
  
  hh(1)=subplot(2,1,1); cla;
  voc_p = waveform(max(1,samp_s-buffer_s):min(samp_e+buffer_e,length(waveform)));
  plot((1:length(voc_p))./Fs,voc_p);
  hold on;
  plot([buffer_s buffer_s]./Fs,[min(voc_p) max(voc_p)],'r')
  plot([buffer_s+samp_e-samp_s buffer_s+samp_e-samp_s]./Fs,...
    [min(voc_p) max(voc_p)],'r')
  axis tight;
  aa=axis;
  axis([aa(1) aa(1)+max_dur+(buffer_s+buffer_e)./Fs aa(3:4)]);
  hold off;
  
  hh(2)=subplot(2,1,2); cla;
  [~,F,T,P] = spectrogram(voc_p,128,120,512,Fs);
  imagesc(T,F,10*log10(P)); set(gca,'YDir','normal');
  set(gca,'clim',[-95 -30]);
  %         colormap jet
  %         colorbar
  hold on;
  axis tight;
  aaa=axis;
  plot([buffer_s buffer_s]./Fs,[0 Fs/2],'r')
  plot([buffer_s+samp_e-samp_s buffer_s+samp_e-samp_s]./Fs,...
    [0 Fs/2],'r')
  hold off;
  
  linkaxes(hh,'x');
  
  reply = getkey;
  while ischar(reply) || isempty(~find(reply==[13 45 27 28 29]))
    disp('neither return or minus were pressed, ESC to quit');
    reply = getkey;
  end
  switch reply
    case {13, 29} %return or right arrow, keep voc
      new_duration_data(vv_indx(vv),2:3)=...
        dur_data(vv_indx(vv),2:3);
      vv=vv+1;
      voc_status(hh,buffer_s/Fs,'OK','g',.05)
    case 45 % delete voc
      new_duration_data(vv_indx(vv),2:3)=nan;
      vv=vv+1;
      voc_status(hh,buffer_s/Fs,'X','r',.15)
    case 27 %ESC
      disp(['On voc: ' num2str(vv_indx(vv))]);
      return;
    case 28 %going backwards
      vv=vv-1;
  end
  if vv < 1
    vv=1;
  elseif vv > length(vv_indx)
    fprintf('<strong>Do you want to save?</strong>\n')
    disp('ESC to cancel, BACK to go back, any other key to continue')
    reply = getkey;
    if isequal(reply, 28)
      vv=vv-2;
    elseif ~isequal(reply, 27)
      saving=1;
      return;
    end
  end
end
