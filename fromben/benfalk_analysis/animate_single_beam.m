clear;
close all;

mic_data_dir='..\mic_data\';

% mic_proc_dir='..\proc_output\';
% trials=dir([mic_proc_dir 'rousettus_20150825*.mat']);

mic_proc_dir='..\proc_output_eptesicus_new\';
trials=dir([mic_proc_dir 'eptesicus_20150824_*_mic_data_bp_proc.mat']);
checked=1;

% interp_method='rb_natural';  % natural neighbor interpolation
interp_method='rb_rbf';  % radial basis function interpolation'
% gaussfitfcn=1; %1 to use fit function of matlab
diag=0;

save_movie=1;
vid_frate=5;

surface_plot=0;

% figure(1), clf, set(gcf,'pos',[10   145   920   480],'color','white')
figure(1), clf, set(gcf,'pos',[10   145   640   480],'color','white')
% colormap hot
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
  
  calls_with_track = bpp.mic_data.call_idx_w_track;
  call_track_locs = round([bpp.mic_data.call.call_start_idx]...
    /bpp.mic_data.fs*bpp.track.fs);
  
  goodch=1:length(bpp.mic_loc(:,1));
  mic_num = 1:bpp.mic_data.num_ch_in_file;
    
  chk_indx=find(bpp.proc.chk_good_call);
  for call_indx=chk_indx'
    if save_movie
      S=strsplit(trials(tt).name,'_');
      trl_indic = strjoin(S(1:4),'_');
      
      v=VideoWriter(['F:\' trl_indic '_call_' num2str(call_indx) '_bp_across_freq_no_surf.avi'],...
        'uncompressed avi');
      v.FrameRate=vid_frate;
      open(v);
    end
    
    ch_good_loc = ~isnan(bpp.mic_loc(:,1))';  % Check for mics without location data
    mic_to_bat_angle = squeeze(bpp.proc.mic_to_bat_angle(call_indx,:,:));
    
    all_freq=unique([bpp.proc.call_freq_vec{call_indx,:}])/1e3;
    for freq_desired=all_freq(all_freq>25 & all_freq< 80)
      call_dB = nan(1,bpp.mic_data.num_ch_in_file);
      for iM=mic_num
        freq = bpp.proc.call_freq_vec{call_indx,iM};
        [~,fidx] = min(abs(freq-freq_desired*1e3));
        call_dB(iM) = bpp.proc.call_psd_dB_comp_re20uPa_withbp{call_indx,iM}(fidx);
      end
      
      ch_ex_sig = find(isnan(call_dB)); % low quality channel from call extraction function
      
      angle_notnanidx = ~ismember(mic_num,ch_ex_sig) & ch_good_loc;
      
      az = mic_to_bat_angle(angle_notnanidx,1);
      el = mic_to_bat_angle(angle_notnanidx,2);
      
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
      vq(~in_smpl_poly) = NaN;
      vq_norm(~in_smpl_poly) = NaN;
      
      [~,II] = max(vq_norm(:));
      [I,J]=ind2sub(size(vq_norm),II);
      
      elqm = elq;
      elqm(isnan(vq_norm)) = NaN;
      azqm = azq;
      azqm(isnan(vq_norm)) = NaN;
      
      vq_norm=vq_norm-max(vq_norm(:));
      
      %plotting the entire beamshape
      cla;
      if surface_plot
        h=surfc_set_zpos(elqm/pi*180,azqm/pi*180,vq_norm);
        h(1).LineStyle='none';
        h(1).FaceLighting='gouraud';
        h(1).FaceAlpha=.9;
        camlight right

        h(2).LineColor = 'k';
        h(2).LevelStep=3;
        view(-80,60)
        axis([-90 90 -90 90 -30 0]);
        xlabel('el')
        ylabel('az')
      else
        contourf(azqm/pi*180,elqm/pi*180,vq_norm,'LevelStep',3,'LineStyle','none')
        view(-90,90)
        axis([-90 90 -90 90]);
        xlabel('az')
        ylabel('el')
      end
      
      caxis([-30 0]);
      title([num2str(freq_desired,'%2.0f') ' kHz']);
      colorbar;
      
      if save_movie
        %  set(gca, 'gridcolor','k','linewidth',1.5,'GridAlpha',.25)
        writeVideo(v,getframe(gcf));
      else
        drawnow
      end
      
      
      if diag
        figure(2), clf, set(gcf,'pos',[168   696   560   420])
        axesm eckert4;
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        axis off

        geoshow(elqm/pi*180,azqm/pi*180,vq_norm,'displaytype','texturemap');
        scatterm(el/pi*180,az/pi*180,[],[0 0 0],'+')
        textm(el/pi*180,az/pi*180,num2str(mic_num(angle_notnanidx)'),...
          'horizontalalignment','center');
        cc = [-30 0];
        caxis(cc);
        colorbar('southoutside');
        tightmap

        vq_norm_min = min(vq_norm(:));
        contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
        cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');
        figure(1);
      end
    end
    
    if save_movie
      close(v);
      fname=[v.Path v.Filename];
      compress_video(fname,1,1,'..\animate_beam_dirs\animate_across_freq\');
    end
    
  end
end

