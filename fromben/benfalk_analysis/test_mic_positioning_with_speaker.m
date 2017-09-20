clear;

%fake data:
room_sz_x=[-2 2];
room_sz_y=[-2 2];
room_sz_z=[0 2];

speed_of_sound=344;

addpath('F:\small_space_beampattern\analysis')

sp_positions=4:2:30;
for ss=sp_positions
  
  trials=50;
  for tt=1:trials
    
    XYZc=[room_sz_x(1) + range(room_sz_x)*rand(ss,1) ...
      room_sz_y(1) + range(room_sz_y)*rand(ss,1)...
      room_sz_z(1) + range(room_sz_z)*rand(ss,1)];
    
    gen_mic_pos=[room_sz_x(1) + range(room_sz_x)*rand(1,1) ...
      room_sz_y(1) + range(room_sz_y)*rand(1,1)...
      room_sz_z(1) + range(room_sz_z)*rand(1,1)];
    
    tof_ch=distance(XYZc,gen_mic_pos)/(speed_of_sound-1);
    
    %corrupt with noise
    std_dev_noise=.05/(speed_of_sound-1); %5 cm, 1 m/s off in speed of sound
    tof_ch=tof_ch + std_dev_noise*randn(ss,1);
    
    
    %euclidean distance formula
    
    %solve this equation for multiple speaker positions
    % speaker_pos_1 = x1,y1,z1
    % speaker_pos_2 = x2,y2,z2
    
    %  (x1-xm)^2+(y1-ym)^2+(z1-zm)^2=d1m^2
    %- (x2-xm)^2+(y2-ym)^2+(z2-zm)^2=d2m^2
    
    mic_pos=nan(1,3);
    for ch=1%:size(sig,2)-1
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
    mic_infer_loc = zeros(size(tof_ch,2),3);
    for iCH=1:size(tof_ch,2)
      notnanidx = ~isnan(tof_ch(:,ch));
      if sum(notnanidx)<4
        mic_infer_loc(iCH,:) = nan(3,1);
      else
        [mic_infer_loc(iCH,:),resnorm] = lsqnonlin(@(x) myfun_p(x,XYZc(notnanidx,:),tof_ch(notnanidx,iCH),c),x0,lb,ub);
      end
    end
    
    
    D1(tt)=distance(mic_pos,gen_mic_pos);
    D2(tt)=distance(mic_infer_loc,gen_mic_pos);
  end
  
  D1sp(ss)=mean(D1);
  D1spstd(ss)=std(D1);
  
  D2sp(ss)=mean(D2);
  D2spstd(ss)=std(D2);
  
end

D1sp(D1sp==0)=[];
D2sp(D2sp==0)=[];
D1spstd(D1spstd==0)=[];
D2spstd(D2spstd==0)=[];



figure(1), clf;
errorbar(sp_positions,D1sp,D1spstd./sqrt(trials))
hold on
errorbar(sp_positions,D2sp,D2spstd./sqrt(trials))
ylabel('Error')
xlabel('# of speakers')
legend({'System of linear eq.','Nonlinear least square'})
a=axis;
axis([a(1:2) 0 .5])


% figure(1), clf;
% histogram(D1)
% hold on;
% histogram(D2)