%allows freq_desired to be passed into script as a variable
a=[];
vars=who;
clear_vars=setdiff(vars,{'freq_desired','save_movie','bat_type'});
clear(clear_vars{:});

close all;

mic_data_dir='..\mic_data\';

if exist('bat_type','var') && strcmp(bat_type,'rousettus')
  mic_proc_dir='..\proc_output\';
  trials=dir([mic_proc_dir 'rousettus_20150825*.mat']);
else
  mic_proc_dir='..\proc_output_eptesicus_new\';
  trials=dir([mic_proc_dir 'eptesicus_20150824_*_mic_data_bp_proc.mat']);
  checked=1;
end

d2=0; %2d data only

% interp_method='rb_natural';  % natural neighbor interpolation
interp_method='rb_rbf';  % radial basis function interpolation'
% gaussfitfcn=1; %1 to use fit function of matlab
diag=0;

if ~exist('freq_desired','var')
  freq_desired=35; %khz
end

if ~exist('save_movie','var')
  save_movie=0;
end
vid_frate=12;
frames_limit_hard=1;

for tt=1:length(trials)  
  bpp=load([mic_proc_dir trials(tt).name]);
  if d2
    bat=bpp.track.track_smooth(:,[1 2]);
  else
    bat=bpp.track.track_smooth;
  end
  
  if exist('checked','var') && checked
    if exist([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat'],'file')
      bpp_checked=load([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat']);
      bpp.proc.chk_good_call = bpp_checked.proc.chk_good_call;
      bpp.proc.ch_ex=bpp_checked.proc.ch_ex;
    else
      continue
    end
  end
  
  if save_movie
    S=strsplit(trials(tt).name,'_');
    trl_indic = strjoin(S(1:4),'_');
    
    out_fname=[trl_indic '_beampattern_' num2str(freq_desired)];
    v=VideoWriter(['F:\' out_fname],...
      'Uncompressed AVI');
    v.FrameRate=vid_frate;
%     v.Quality=95;
    open(v);
  end
  
  calls_with_track = bpp.mic_data.call_idx_w_track;
  call_track_locs = round([bpp.mic_data.call.call_start_idx]...
    /bpp.mic_data.fs*bpp.track.fs);
  
  goodch=1:length(bpp.mic_loc(:,1));
  mic_num = 1:bpp.mic_data.num_ch_in_file;
  
  if frames_limit_hard
    frames = call_track_locs( calls_with_track( find(bpp.proc.chk_good_call,1) ) ) :...
      call_track_locs( calls_with_track( find(bpp.proc.chk_good_call,1,'last') ) );
    
    %limit to speed > 1.5 m/s
    speed = calc_speed(bat,bpp.track.fs,0);
    speed = [nan; speed];
    good_speed = find(speed > 1.5);
    good_speed (good_speed < frames(1))=[];
    good_speed (good_speed > frames(end))=[];
    
    frames = max(good_speed(1),frames(1)) : min(good_speed(end),frames(end));
  else
    frames=find(isfinite(bat(:,1)));
    frames=max(frames(1),...
      call_track_locs(calls_with_track( find(bpp.proc.chk_good_call,1) ) )-10):...
      min(call_track_locs(calls_with_track(find(bpp.proc.chk_good_call,1,'last')))+10,...
      frames(end));
  end
  figure(1); clf; set(gcf,'pos',[10 40 800 636],'color','w')
  axesm eckert4;
  framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
  gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
  axis off
  tightmap
  cb=colorbar('location','southoutside');
  set(cb,'visible','off')
  
  for fr=frames(1):frames(end)
    call_present=ismember(call_track_locs(calls_with_track),fr);
    call_indx=find(call_present);
    if ~isempty(call_indx) && ismember(fr,frames)    
      call_dB = nan(1,bpp.mic_data.num_ch_in_file);
      for iM=mic_num
        freq = bpp.proc.call_freq_vec{call_indx,iM};
        [~,fidx] = min(abs(freq-freq_desired*1e3));
        call_dB(iM) = bpp.proc.call_psd_dB_comp_re20uPa_withbp{call_indx,iM}(fidx);
      end
      
      ch_ex_sig = find(isnan(call_dB)); % low quality channel from call extraction function
      ch_good_loc = ~isnan(bpp.mic_loc(:,1))';  % Check for mics without location data
      
      angle_notnanidx = ~ismember(mic_num,ch_ex_sig) & ch_good_loc;
      
      mic_to_bat_angle = squeeze(bpp.proc.mic_to_bat_angle(call_indx,:,:));
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
      
      %plotting the entire beamshape
      cla
      axesm eckert4;
      framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
      gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
      axis off
      geoshow(elqm/pi*180,azqm/pi*180,vq_norm,'displaytype','texturemap');
      contourm(elq/pi*180,azq/pi*180,vq_norm,-3,'w','linewidth',2);
      scatterm(el/pi*180,az/pi*180,[],[0 0 0],'+')
%       textm(el/pi*180,az/pi*180,num2str(mic_num(angle_notnanidx)'),'horizontalalignment','center');
      cc = [-30 0];
      caxis(cc);
      colorbar('southoutside');
%       tightmap
%       setm(gca,'ParallelLabel','on','MeridianLabel','on')

      vq_norm_min = min(vq_norm(:));
      contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
      cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');


      if diag
      end
      
      if ~save_movie
        drawnow
        pause(.15)
      end
    else
      
    end
    
    if save_movie
%       set(gca, 'gridcolor','k','linewidth',1.5,'GridAlpha',.25)
      writeVideo(v,getframe(gcf));
    else
      drawnow
    end
  end
  
  if save_movie
    close(v);
    
    [~,~]=system(['C:\video_tools\ffmpeg_64\bin\ffmpeg -y ' ...
      '-i F:\' out_fname '.avi '...
      '-c:v libx264 -crf 20 -pix_fmt yuv420p ' ...
      '..\animate_beam_dirs\' out_fname '.mp4']);
    delete(['F:\' out_fname '.avi'])
  end
end