function extract_dur2file(fn)

DIAG=0;

% pname='..\mic_recordings\20150911\eptesicus_LB53_matfile\';
% fn='eptesicus_LB53_1_mic_data_detect.mat';

if ( nargin < 1 || isempty(fn) )
  if ispref('small_space_beam') && ispref('small_space_beam','micrecpath')
    pname=getpref('small_space_beam','micrecpath');
  else
    pname='';
  end
  [fname,pname]=uigetfile('*.mat',[],pname);
  if isequal(fname,0)
    return;
  end
  fn=[pname,fname];
  setpref('small_space_beam','micrecpath',pname)
end

load(fn)

ch=mode([call.channel_marked]);
Y=sig(:,ch);
pretrig_t=length(Y)/fs;

voc_s=[call([call.channel_marked]==ch).locs]';
voc_t=-pretrig_t+voc_s./fs;

[onsets,offsets,voc_t]=extract_dur(Y,fs,voc_t,-pretrig_t,0,pretrig_t,[],0,DIAG);
dur_data=[voc_t, onsets, offsets];
new_dur_data=dur_data;
[saving,new_dur_data]=mark_good_dur(dur_data,Y,fs,pretrig_t,new_dur_data);
checked=1;
manual_marked=0;
timestamp=now;

if saving
  [pname,fname] = fileparts(fn);
  
  save([pname '\' fname '_processed_duration.mat'],...
    'new_dur_data','pretrig_t','ch','fs','checked','manual_marked','timestamp',...
    'fname','pname')
end



