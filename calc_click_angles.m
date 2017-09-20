%calc click angles
clear; close all;

mic_proc_dir='..\proc_output\';
trials=dir([mic_proc_dir 'rousettus_20150825*.mat']);
d2=0; %2d data only
% interp_method='rb_natural';  % natural neighbor interpolation
interp_method='rb_rbf';  % radial basis function interpolation'
% gaussfitfcn=1; %1 to use fit function of matlab
use_interp_beamshape=0;
DIAG=0;

freq_desired=35; %khz

angle_diff=[];

for tt=1:length(trials)
  bpp=load([mic_proc_dir trials(tt).name]);
  if d2
    bat=bpp.track.track_smooth(:,[1 2]);
  else
    bat=bpp.track.track_smooth;
  end
  
  calls_with_track = bpp.mic_data.call_idx_w_track;
  call_track_locs = round([bpp.mic_data.call.call_start_idx]...
    /bpp.mic_data.fs*bpp.track.fs);
  call_w_track_locs=call_track_locs(calls_with_track);
  good_calls = logical(bpp.proc.chk_good_call);
  
  goodch=1:length(bpp.mic_loc(:,1));
  mic_num = 1:bpp.mic_data.num_ch_in_file;
  
  beam_dirs=nan(size(find(good_calls)',1),1);
  for call_indx=find(good_calls)'
    fr=call_w_track_locs(call_indx);
    
    call_dB = nan(1,bpp.mic_data.num_ch_in_file);
    for iM=mic_num
      freq = bpp.proc.call_freq_vec{call_indx,iM};
      [~,fidx] = min(abs(freq-freq_desired*1e3));
      call_dB(iM) = bpp.proc.call_psd_dB_comp_re20uPa_withbp{call_indx,iM}(fidx);
    end
    
    ch_ex_sig = find(isnan(call_dB)); % low quality channel from call extraction function
    ch_good_loc = ~isnan(bpp.mic_loc(:,1))';  % Check for mics without location data
    
    angle_notnanidx = ~ismember(mic_num,ch_ex_sig) & ch_good_loc;
    
    %create my own angles for bat to mic
    mic_vec=bpp.mic_loc(goodch,:)-...
      repmat(bat(fr,:),size(bpp.mic_loc(goodch,:),1),1);
    [az,el] =...
      cart2sph(mic_vec(:,1),mic_vec(:,2),mic_vec(:,3));
    
    %       mic_to_bat_angle = squeeze(bpp.proc.mic_to_bat_angle(call_indx,:,:));
    %       az = mic_to_bat_angle(angle_notnanidx,1);
    %       el = mic_to_bat_angle(angle_notnanidx,2);
    
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
    
%     head_dir_xyz=bpp.head_aim.head_aim_int(...
%       bpp.track.call_loc_idx_on_track_interp(call_indx),:);
%     [h_az,h_el]=cart2sph(head_dir_xyz(1),head_dir_xyz(2),head_dir_xyz(3));
    
%     [bx,by,bz]=sph2cart(azqm(I,J),elqm(I,J),.2);
    
    
    beam_dirs(call_indx) = azqm(I,J);
    
    if DIAG
      %plotting the beam direction
      if bpp.proc.chk_good_call(call_indx)
        beam_col=[55 126 184]./255;
      else
        beam_col=[228 26 28]./255;
      end
      plot3([bat(fr,1) bat(fr,1)+bx],...
        [bat(fr,2) bat(fr,2)+ by],...
        [bat(fr,3) bat(fr,3)+ bz],...
        '-','linewidth',2,'color',beam_col);
    end
  end
  
  callt=[bpp.mic_data.call.call_start_idx]./bpp.mic_data.fs;
  callt_good=callt(calls_with_track(good_calls));
  PI=diff(callt_good);
  beam_dirs_good=beam_dirs(good_calls);
  
  pairs = PI < .03;
  for pp=find(pairs)
    pp1=pp;
    pp2=pp+1;
    angle_diff(end+1) = abs(beam_dirs_good(pp1)-beam_dirs_good(pp2));
  end
end
