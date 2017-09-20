%calc click angles
clear; 

mic_proc_dir='..\proc_output\';
species_date='rousettus_20150825';
trials=dir([mic_proc_dir species_date '*.mat']);
d2=0; %2d data only
% interp_method='rb_natural';  % natural neighbor interpolation
interp_method='rb_rbf';  % radial basis function interpolation'
% gaussfitfcn=1; %1 to use fit function of matlab
use_interp_beamshape=0;
DIAG=0;

freq_desired=35; %khz

[angle_diff,beam_dirs_az,beam_dirs_el,...
  bat_dist2mics,voc_t,callnums,bat_pos,beam_I]=deal(cell(length(trials),1));

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
  
  beam_dirs_az{tt}=nan(size(good_calls,1),1);
  beam_dirs_el{tt}=nan(size(good_calls,1),1);
  bat_dist2mics{tt}=nan(size(good_calls,1),1);
  voc_t{tt}=nan(size(good_calls,1),1);
  callnums{tt}=nan(size(good_calls,1),1);
  bat_pos{tt}=nan(size(good_calls,1),3);
  beam_I{tt}=nan(size(good_calls,1),1);
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
    
    [azq,elq] = meshgrid(-pi/2:pi/180:pi/2,-pi/2:pi/180:pi/2);
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
    
    beam_dirs_az{tt}(call_indx) = azqm(I,J);
    beam_dirs_el{tt}(call_indx) = elqm(I,J);
    bat_dist2mics{tt}(call_indx)=min(distance(bpp.mic_loc(goodch,:),bat(fr,:)));
    voc_t{tt}(call_indx)=bpp.mic_data.call(call_indx).call_start_idx./bpp.mic_data.fs;
    callnums{tt}(call_indx)=call_indx;
    bat_pos{tt}(call_indx,:)=bat(fr,:);
    beam_I{tt}(call_indx)=max(vq(:));
    
    if DIAG
      %plotting the beam direction
      figure(1); clf;
      if bpp.proc.chk_good_call(call_indx)
        beam_col=[55 126 184]./255;
      else
        beam_col=[228 26 28]./255;
      end
      [bx,by,bz]=sph2cart(beam_dirs_az{tt}(call_indx),beam_dirs_el{tt}(call_indx),.5);
      plot3([bat(fr,1) bat(fr,1)+bx],...
        [bat(fr,2) bat(fr,2)+ by],...
        [bat(fr,3) bat(fr,3)+ bz],...
        '-','linewidth',2,'color',beam_col);
      hold on;
      plot3(bat(fr,1),bat(fr,2),bat(fr,3),'or');
      scatter3(bpp.mic_loc(:,1),bpp.mic_loc(:,2),bpp.mic_loc(:,3))
      axis equal; grid on; view(2); pause(.1)
    end
  end
  
  PI=diff(voc_t{tt})';
  pairs = PI < .03;
  for pp=find(pairs)
    pp1=pp;
    pp2=pp+1;
    
    [p1x,p1y,p1z]=sph2cart(beam_dirs_az{tt}(pp1),beam_dirs_el{tt}(pp1),1);
    [p2x,p2y,p2z]=sph2cart(beam_dirs_az{tt}(pp2),beam_dirs_el{tt}(pp2),1);
    
    p1vec=[p1x,p1y,p1z];
    p2vec=[p2x,p2y,p2z];
    
    angle_diff{tt}(end+1)=atan2(norm(cross(p1vec,p2vec)),dot(p1vec,p2vec));
  end
end

save(['extracted_beam_dirs_' species_date '.mat'],'angle_diff','beam_dirs_az','beam_dirs_el',...
  'bat_dist2mics','voc_t','callnums','bat_pos','trials','beam_I')