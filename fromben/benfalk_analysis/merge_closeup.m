clear;

aud_dir='..\mic_recordings\20151218_rousettus_matfile\';
miro_dir='..\miro\20151218\';
merged_dir='..\closeup_merged\';

afiles=dir([aud_dir '*.mat']);
afnames={afiles.name};
CC=cellfun(@(c) strsplit(c,'_'),afnames,'uniformoutput',0);
[~,isort]=sort( cell2mat(cellfun(@(c) str2double(c{4}),CC,'uniformoutput',0)') );
afnames=afnames(isort)';
upsamp_factor=5;

% mfiles=dir([miro_dir '*.mp4']);
cams={'c_14785','c_14786','c_16645'};
cam_files=cellfun(@(c) dir([miro_dir c '*.mp4']),cams,'uniformoutput',0);

load('miro_mic_20151218.mat'); %miro_mic_20151218
frm_rts=[];
frm_rts(1:length(miro_mic_20151218))=500;
frm_rts(11)=2e3;

for ff=15%1:size(miro_mic_20151218,1)
  mmidx=miro_mic_20151218(ff,1);
  aaidx=miro_mic_20151218(ff,2);
  
  if isnan(aaidx) || isnan(mmidx)
    continue
  end
  
  load([aud_dir afnames{aaidx}])
  [MM,idx]=max(max(abs(sig))./median(abs(sig)));
  
  aud_sig = sig(:,idx)./max(abs(sig(:,idx)));
  resamp_aud = resample(aud_sig,upsamp_factor,1);
  fs=upsamp_factor*fs;
  
  audiowrite([merged_dir afnames{aaidx}(1:end-4) '.flac'],...
    resamp_aud./max(abs(resamp_aud)),fs/(frm_rts(ff)/25));
  
  for mc=1:3
    camfile=cam_files{mc}(mmidx).name;
    [status,cmdout]=system(['C:\video_tools\ffmpeg_64\bin\ffmpeg.exe -y'...
      ' -i ' miro_dir camfile ...
      ' -i ' merged_dir afnames{aaidx}(1:end-4) '.flac' ...
      ' -c copy'...
      ' ' merged_dir camfile(1:end-4) '_merged.mkv']);
  end
  
end