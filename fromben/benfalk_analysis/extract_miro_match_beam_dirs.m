clear; close all;

% mic_proc_dir='..\proc_output\';
% trials=dir([mic_proc_dir 'rousettus_20150825*.mat']);

mic_proc_dir='..\proc_output_eptesicus_new\';
trials=dir([mic_proc_dir 'eptesicus_20150824_*_mic_data_bp_proc.mat']);
checked=1;

out_dir = '..\animate_beam_dirs\';

vid_frate=12;
frames_limit_hard=1;

for tt=1:length(trials)
  bpp=load([mic_proc_dir trials(tt).name]);
  
  if exist('checked','var') && checked 
    if exist([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat'],'file')
      bpp_checked=load([mic_proc_dir trials(tt).name(1:end-4) '_checked.mat']);
      bpp.proc.chk_good_call = bpp_checked.proc.chk_good_call;
      bpp.proc.ch_ex=bpp_checked.proc.ch_ex;
    else
      continue
    end
  end
  
  bat=bpp.track.track_smooth;
  
  calls_with_track = bpp.mic_data.call_idx_w_track;
  call_track_locs = round([bpp.mic_data.call.call_start_idx]...
    /bpp.mic_data.fs*bpp.track.fs);
  
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
      call_track_locs(calls_with_track(find(bpp.proc.chk_good_call,1)))-15):...
      min(call_track_locs(calls_with_track(find(bpp.proc.chk_good_call,1,'last')))+15,...
      frames(end));
  end
  
  startf=frames(1);
  endf=frames(end);
  
  mfnames=get_mfname_from_afname(trials(tt).name);
  cnames_split=cellfun(@(c) strsplit(c,{'_','\\'}),mfnames,'uniformoutput',0);
  tname=strsplit(trials(tt).name,'_');
  
  for mm=1:length(mfnames)
    [status,cmdout]=system(['C:\video_tools\ffmpeg_64\bin\ffmpeg.exe -y '...
      ' -r "' num2str(vid_frate) '"'...
      ' -i ' mfnames{mm} ...
      ' -ss ' num2str((startf-1)/vid_frate) ... 
      ' -t ' num2str((endf-startf+1)/vid_frate) ' -c:v libx264 -crf 20 '...
      ' -pix_fmt yuv420p '...
      ' ' out_dir strjoin(tname(1:4),'_') '_c' cnames_split{mm}{5} '.mp4']);
  end
end