clear;

diag=0;

base_dir='F:\small_space_beampattern\analysis\';

mic_data_path='..\mic_data\';
datadate='20150826';
mic_recs=dir([base_dir mic_data_path 'speaker*' datadate '*calib*.mat']);

num_ch=34;

time_btw_pulses=.5; %sec
pulse_width=.003;
max_tof = .01;

speed_of_sound = 344; %50% humidity, 20 C, approximated

load([base_dir 'speaker_pos_20150826_circle_fit.mat']);


tof_ch=nan(length(mic_recs),num_ch);
for k=1:length(mic_recs)
  
  load([base_dir mic_data_path mic_recs(k).name])
  
  for ch=1:size(sig,2)-1
    search_sig=sig(:,end);
    ch_sig = sig(:,ch);
    
    search_sig_square=search_sig.^2;
    [pks,locs]=findpeaks(search_sig_square,fs,'minpeakheight',max(search_sig_square)/2,...
      'minpeakdistance',time_btw_pulses-.05); %minpeakdistance in time (s)
    
    if diag
      figure(1); clf;
      t=((1:length(search_sig))-1)./fs;
      plot(t,search_sig)
      hold on;
      plot(t,ch_sig);
      title(['Channel ' num2str(ch)])
    end
    
    tof=nan(length(locs),1);
    for ll=1:length(locs)
      sig_win=round([locs(ll)-max_tof,locs(ll)+max_tof]*fs);
      sig_ex = search_sig(sig_win(1):sig_win(2));
      
      data_win=round([locs(ll)-max_tof,locs(ll)+max_tof]*fs);
      data_ex = ch_sig(data_win(1):data_win(2));
      
      [xc,lags] = xcorr(sig_ex,data_ex);
      
      xc_envel=abs(hilbert(xc));
      [~,indx]=max(xc_envel);
      
      tof(ll)=abs(lags(indx))/fs;
      
      if diag
        figure(2); set(gcf, 'pos', [50 50 450 950]); clf;
        axh1=subaxis(3,1,1,'m',.075,'mt',.04,'ml',.12);
        plot((1:length(sig_ex))/fs*1e3,sig_ex);
        hold on;
        plot((1:length(sig_ex))/fs*1e3,data_ex);
        plot(max_tof*1e3,max(sig_ex)*.9,'*r')
        plot(max_tof*1e3 - lags(indx)/fs*1e3,max(sig_ex)*.9,'*g')
        title(['Channel ' num2str(ch)])
        axis tight;
        
        axh2=subaxis(3,1,2,'m',.075,'mt',.04,'ml',.12);
        sig_comb=sig_ex+data_ex;
        spectrogram(sig_comb,256,250,256,fs,'yaxis');
        colorbar off
        hold on;
        plot(repmat(max_tof*1e3,2,1),[0 fs/2/1e3],'r')
        plot( repmat(( max_tof - lags(indx)/fs )*1e3,2,1),[0 fs/2/1e3],'g')
        
        linkaxes([axh1 axh2],'x')
        
        subaxis(3,1,3,'m',.075,'mt',.04,'ml',.12)
        plot(xc_envel)
        hold on;
        plot(abs(xc))
        plot(indx,max(xc_envel),'*r')
        a=axis;
        axis([indx-pulse_width*fs/2 indx+pulse_width*fs/2 a(3:4)])
        
        drawnow;
        pause(.01);
      end
    end
    
    tof_ch(k,ch)=nanmean(tof);
  end
  
end

%euclidean distance formula

%solve this equation for multiple speaker positions
% speaker_pos_1 = x1,y1,z1
% speaker_pos_2 = x2,y2,z2

%  (x1-xm)^2+(y1-ym)^2+(z1-zm)^2=d1m^2
%- (x2-xm)^2+(y2-ym)^2+(z2-zm)^2=d2m^2

mic_pos=nan(size(sig,2)-1,3);
for ch=1:size(sig,2)-1
  s_x=XYZc(:,1);
  s_y=XYZc(:,2);
  s_z=XYZc(:,3);
  
  d1m=tof_ch(1,ch)*speed_of_sound;
  x1=s_x(1);
  y1=s_y(1);
  z1=s_z(1);
  for m=1:length(XYZc)-1
    sp_pos=m+1;
    d2m=tof_ch(sp_pos,ch)*speed_of_sound;
    x2=s_x(sp_pos);
    y2=s_y(sp_pos);
    z2=s_z(sp_pos);
    
    C(m)=d1m^2-d2m^2+x2^2-x1^2+y2^2-y1^2+z2^2-z1^2;
    B(m,:)=[2*(-x1+x2), 2*(-y1+y2), 2*(-z1+z2)];
  end
  
	mic_pos(ch,:)=(B\C')';
end


% Infer mic locations
%from wu-jung
x0 = [1,1,1];
ub = [5,5,5]; lb = [-5,-5,-5];
c = 344;
mic_infer_loc = zeros(length(tof_ch),3);
for iCH=1:length(tof_ch)
  notnanidx = ~isnan(tof_ch(:,ch));
  if sum(notnanidx)<4
    mic_infer_loc(iCH,:) = nan(3,1);
  else
    [mic_infer_loc(iCH,:),resnorm] = lsqnonlin(@(x) myfun_p(x,XYZc(notnanidx,:),tof_ch(notnanidx,iCH),c),x0,lb,ub);
  end
end


if diag
  figure(4), clf;
  h=plot3(mic_pos(:,3),mic_pos(:,1),mic_pos(:,2),'+');
  text(mic_pos(:,3),mic_pos(:,1),mic_pos(:,2),num2str( (1:length(mic_pos) )' ),...
    'color',h.Color)
  axis equal, grid on; hold on;
  
  h=plot3(mic_infer_loc(:,3),mic_infer_loc(:,1),mic_infer_loc(:,2),'x');
  text(mic_infer_loc(:,3),mic_infer_loc(:,1),mic_infer_loc(:,2),num2str( (1:length(mic_infer_loc) )' ),...
    'color',h.Color)
  
%   MP_speaker_old=load([base_dir 'mic_pos_via_speaker_20150826.mat']);
%   hold on;
%   h=plot3(MP_speaker_old.mic_infer_loc(:,3),MP_speaker_old.mic_infer_loc(:,1),...
%     MP_speaker_old.mic_infer_loc(:,2),'o');
%   text(MP_speaker_old.mic_infer_loc(:,3),MP_speaker_old.mic_infer_loc(:,1),...
%     MP_speaker_old.mic_infer_loc(:,2),...
%     num2str( (1:length(MP_speaker_old.mic_infer_loc) )' ),'color',h.Color)
  
  MP=load([base_dir '..\mic_pos\mic_pos_20150826_S1_12.mat']);
  hold on;
  tip_indx=[(1:3:32*3)+2 32*3+2 32*3+4];
  tips=MP.mic_pos(tip_indx,:);
  h=plot3(tips(:,3),tips(:,1),tips(:,2),'s');
  text(tips(:,3),tips(:,1),tips(:,2),num2str((1:length(tips(:,1)))'),'color',h.Color)
  
  legend({'Ben (system of linear eq.)','Wu-Jung (nonlinear)','Vicon'},'location','southoutside','orientation','horizontal')
  
  
  addpath('F:\small_space_beampattern\analysis')
  [mic_err_mic_pos, mic_err_mic_infer_loc]=deal(nan(length(mic_pos),1));
  for mm=1:length(mic_pos)
    DD=distance(mic_pos(mm,:),tips);
    mic_err_mic_pos(mm)=min(DD);
    
    DD=distance(mic_infer_loc(mm,:),tips);
    mic_err_mic_infer_loc(mm)=min(DD);
  end
  mean(mic_err_mic_pos)
  std(mic_err_mic_pos)
  mean(mic_err_mic_infer_loc)
  std(mic_err_mic_infer_loc)
end

