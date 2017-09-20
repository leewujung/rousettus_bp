% 2015 12 03  Merge data from all clicks
% 2015 12 12  Use 35 kHz rotation for other frequencies
% 2015 12 22  Adopt the merging click code for finding multi-freq center of
%             the beampatterns
% 2015 12 24  Move the plotting part to 'compile_info_plot.m'
%             This is only for processing
% 2015 12 24  Separate out the shift and rotation portion and save data

clear
warning off

bat='36134';

colordef black

base_path = '..\..\';
save_path = fullfile(base_path,'proc_output_beam_align');
if ~exist(save_path,'dir')
  mkdir(save_path);
end

%bat_proc_path = '.\proc_output_eptesicus_new';
%bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,['eptesicus_20150824_',bat,'*bp_proc.mat']));
%checked=1;

 bat_proc_path = '.\proc_output';
 bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,['rousettus_20150825_',bat,'*bp_proc.mat']));

data_fname = ['all_click_rotated_' bat '.mat'];
plot_opt_all_click = 1;
plot_opt_indiv_click = 0;
save_opt = 1;
freq_rot = 35e3;

addpath('F:\repositories\beampattern_processing\fit_beams\');
addpath('..\')

%requires:
% addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
% addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
% addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);

A.processed_files = bat_proc_file_all;
A.base_path = base_path;
A.data_path = bat_proc_path;
A.save_path = save_path;
A.freq_used_for_rotation = freq_rot;

% Set map projection
map_proj = 'eckert4';
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
[xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
xxx = xxx(1:2);
yyy = yyy(3:4);

A.map_projection = map_proj;
A.mstruct = mstruct;
A.map_plot_xlim = xxx;
A.map_plot_ylim = yyy;


% Process all files
tic
raw_meas.az = cell(length(bat_proc_file_all),1);
raw_meas.el = cell(length(bat_proc_file_all),1);
raw_meas.call_dB = cell(length(bat_proc_file_all),1);
raw_meas.ch_include_idx = cell(length(bat_proc_file_all),1);
raw_meas.click_side = cell(length(bat_proc_file_all),1);
raw_meas.click_side_rep = cell(length(bat_proc_file_all),1);
raw_meas.click_num = cell(length(bat_proc_file_all),1);
rotate_data = cell(length(bat_proc_file_all),1);
shift_tilt_final.az = cell(length(bat_proc_file_all),1);
shift_tilt_final.el = cell(length(bat_proc_file_all),1);
for iB = 1:length(bat_proc_file_all)
  
  bat_proc_file = bat_proc_file_all(iB).name;
  ss = strsplit(bat_proc_file,'_');
  save_fname = strjoin(ss(3:4),'-');
  
  fprintf('File: %s\n',bat_proc_file);
  
  data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
  if exist('checked','var') && checked 
    if exist(fullfile(base_path,bat_proc_path,[bat_proc_file_all(iB).name(1:end-4)...
        '_checked.mat']),'file')
      bpp_checked=load(fullfile(base_path,bat_proc_path,[bat_proc_file_all(iB).name(1:end-4)...
        '_checked.mat']));
      data.proc.chk_good_call = bpp_checked.proc.chk_good_call;
      data.proc.ch_ex=bpp_checked.proc.ch_ex;
    else
      continue
    end
  end
  
  good_call_idx = find(data.proc.chk_good_call);
  
  if plot_opt_all_click
    numrow = ceil(length(good_call_idx)/4);
    fig_clicks = figure;
    set(fig_clicks,'Position',[5 100 1000 140*numrow]);
  end
  
  for iC = good_call_idx'
    
    iC_save = find(iC==good_call_idx);
    fprintf('Call %02d\n',iC);
    
    % Get call info and rotate beampattern using 35 kHz beampattern
    [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_rot,iC);
    [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt,fig_fit] = ...
      shift_rotate_bp(az(ch_include_idx),el(ch_include_idx),call_dB(ch_include_idx),...
      'eckert4',plot_opt_indiv_click);
    
    if plot_opt_indiv_click
      figure(fig_fit);
      suptitle(sprintf('%s, Call #%02d',save_fname,iC));
      saveSameSize(fig_fit,'file',...
        fullfile(save_path,sprintf('indiv_click_file_%s_call%d.png',save_fname,iC)),...
        'format','png','renderer','painters');
      close(fig_fit)
    end
    
    raw_meas.click_num{iB}(iC_save,:) = iC*ones(length(az),1);
    raw_meas.az{iB}(iC_save,:) = az;
    raw_meas.el{iB}(iC_save,:) = el;
    raw_meas.call_dB{iB}(iC_save,:) = call_dB;
    raw_meas.ch_include_idx{iB}(iC_save,:) = ch_include_idx;
    raw_meas.click_side{iB}(iC_save) = click_side;
    raw_meas.click_side_rep{iB}(iC_save,:) = click_side*ones(length(az),1);
    
    shift_tilt_final.az{iB}(iC_save,:) = rot_elpctr_tilt.az;  % final shifted & rotated azimuth
    shift_tilt_final.el{iB}(iC_save,:) = rot_elpctr_tilt.el;  % final shifted & rotated elevation
    
    rotate_data{iB}(iC_save).raw = raw;
    rotate_data{iB}(iC_save).ror_max = rot_max;
    rotate_data{iB}(iC_save).rot_elpctr = rot_elpctr;
    rotate_data{iB}(iC_save).rot_elpctr_tilt = rot_elpctr_tilt;
    
    % Plot rotated bp using ellipse center
    if plot_opt_all_click
      figure(fig_clicks)
      subaxis(numrow,4,find(iC==good_call_idx),'m',.05);
      contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,raw.vq_norm,-3:-3:-39,'fill','on');
      axis equal
      grid on
      xlim([xxx(1) xxx(2)])
      ylim([yyy(1) yyy(2)])
      caxis([-20 0])
      title(sprintf('Call #%02d',iC));
      drawnow
    end
  end
  
  if plot_opt_all_click
    figure(fig_clicks)
    set(gcf,'color','black')
    suptitle([save_fname,' all']);
    export_fig(fullfile(save_path,['all_clicks_',save_fname]),'-png','-nocrop')
    close(fig_clicks)
  end
  
end
warning on
toc

A.raw_meas = raw_meas;
A.rotate_data = rotate_data;
A.shift_tilt_final = shift_tilt_final;

if save_opt==1
  disp('Saving...')
  save(fullfile(save_path,data_fname),'-struct','A');
  disp('Saved')
end
