% 2015 12 03  Merge data from all clicks
% 2015 12 12  Use 35 kHz rotation for other frequencies
% 2015 12 22  Adopt the merging click code for finding multi-freq center of
%             the beampatterns
% 2015 12 24  Move the plotting part to 'compile_info_plot.m'
%             This is only for processing
% 2015 12 24  Separate out the shift and rotation portion and save data
% 2016 04 19  Deal with "out-of-boundary" problem
%             updated shift_rotate_bp.m and bp_fit_ellipse_azel.m
% 2016 05 07  Eptesicus beams
clear
warning off
usrn = getenv('username');

if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

% Params and options
plot_opt_all_click = 1;
plot_opt_indiv_click = 1;
save_opt = 1;

cvec = 0:-3:-39;
freq_rot = 35e3;
it_shift_th = 0.005;
azel_bnd = [300,150];

% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
end

save_root = 'analysis_results_figs';
save_path = fullfile(base_path,save_root,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

bat_proc_path = 'proc_output_eptesicus_new';
bat = {'LB53'};

for iBAT=1:length(bat)
if strcmp(bat{iBAT},'3bat')
   bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'eptesicus_20150825_*'));
else
   bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,['eptesicus_20150825_',bat{iBAT},'*']));
end
save_fname = ['all_click_rotated_',bat{iBAT},'.mat'];

A.base_path = base_path;
A.save_path = save_path;
A.bp_processed_path = bat_proc_path;
A.bp_processed_file_all = bat_proc_file_all;
A.freq_used_for_rotation = freq_rot;
A.iterative_shift_threshold = it_shift_th;
A.shift_bound_azel = azel_bnd;

% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
[xxx,yyy] = mfwdtran(mstruct,[0,0,-90,90],[-180,180,0,0]);
xxx = xxx(1:2);
yyy = yyy(3:4);
xy_lim = [xxx(1) xxx(2) yyy(1) yyy(2)];
    
A.map.map_projection = map_proj;
A.map.mstruct = mstruct;
A.map.map_plot_xlim = xxx;
A.map.map_plot_ylim = yyy;
A.map.map_plot_xylim = xy_lim;

% Process all files
tic
raw_meas.az = cell(length(bat_proc_file_all),1);
raw_meas.el = cell(length(bat_proc_file_all),1);
raw_meas.call_dB = cell(length(bat_proc_file_all),1);
raw_meas.ch_include_idx = cell(length(bat_proc_file_all),1);
raw_meas.click_side = cell(length(bat_proc_file_all),1);
raw_meas.click_side_rep = cell(length(bat_proc_file_all),1);
raw_meas.click_num = cell(length(bat_proc_file_all),1);
raw_meas.good_call_idx = cell(length(bat_proc_file_all),1);
rotate_data = cell(length(bat_proc_file_all),1);
shift_tilt_final.az = cell(length(bat_proc_file_all),1);
shift_tilt_final.el = cell(length(bat_proc_file_all),1);
shift_tilt_final.outofbnd_azel_idx = cell(length(bat_proc_file_all),1);
 for iB = 1:length(bat_proc_file_all)

    bat_proc_file = bat_proc_file_all(iB).name;
    ss = strsplit(bat_proc_file,'_');
    save_fname = strjoin([script_name,ss(3:4)],'_');

    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    good_call_idx = find(data.proc.chk_good_call);
    raw_meas.good_call_idx{iB} = good_call_idx;
    
    if plot_opt_all_click
        numrow = ceil(length(good_call_idx)/4);
        fig_clicks = figure;   % shift/tilt
        set(fig_clicks,'Position',[100 100 1000 140*numrow]);
        fig_clicks_ofbfix = figure;   % shift/tilt and delete out-of-bnd points
        set(fig_clicks_ofbfix,'Position',[100 100 1000 140*numrow]);
    end
    
    for iC = good_call_idx'

        iC_save = find(iC==good_call_idx);
        fprintf('Call %02d\n',iC);
        
        % Get call info and rotate beampattern using 35 kHz beampattern
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_rot,iC);
        [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
            shift_rotate_bp(az(ch_include_idx),el(ch_include_idx),call_dB(ch_include_idx),'eckert4',it_shift_th,azel_bnd);
        
        % Check az, el shift and tilt
        if plot_opt_indiv_click==1
            Vpass.raw = raw;
            Vpass.rot_max = rot_max;
            Vpass.rot_elpctr = rot_elpctr;
            Vpass.rot_elpctr_tilt = rot_elpctr_tilt;
            
            % Movement of [az,el] of each mic
            fig_azel = plot_indiv_click_azel_movement(Vpass);
            title(sprintf('%s, Call #%02d',regexprep(save_fname,'_','\\_'),iC));
            saveSameSize(fig_azel,'file',...
                fullfile(save_path,sprintf('%s_c%02d_azel_chk.png',save_fname,iC)),...
                'format','png','renderer','painters');
            close(fig_azel)

            % Plot the shift/tilt procedure
            fig_fit = plot_indiv_click_rotate(Vpass,cvec,xy_lim);
            suptitle(sprintf('%s, Call #%02d',regexprep(save_fname,'_','\\_'),iC));
            saveSameSize(fig_fit,'file',...
                fullfile(save_path,sprintf('%s_c%02d_contour.png',save_fname,iC)),...
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
        shift_tilt_final.outofbnd_azel_idx{iB}(iC_save,:) = rot_elpctr_tilt.outofbnd_azel_idx;

        rotate_data{iB}(iC_save).raw = raw;
        rotate_data{iB}(iC_save).rot_max = rot_max;
        rotate_data{iB}(iC_save).rot_elpctr = rot_elpctr;
        rotate_data{iB}(iC_save).rot_elpctr_tilt = rot_elpctr_tilt;
        
        % Plot rotated bp using ellipse center
        if plot_opt_all_click
            figure(fig_clicks)
            subplot(numrow,4,find(iC==good_call_idx));
            contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,raw.vq_norm,cvec,'fill','on');
            axis equal
            grid on
            axis(xy_lim)
            title(sprintf('Call #%02d',iC));
            colormap(parula(length(cvec)-1));
            caxis([cvec(end), 0])
            
            figure(fig_clicks_ofbfix)   % delete out-of-bnd points
            subplot(numrow,4,find(iC==good_call_idx));
            vq_norm_fix = raw.vq_norm;
            vq_norm_fix(rot_elpctr_tilt.outofbnd_azel_idx_q) = NaN;
            contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,...
                    vq_norm_fix,cvec,'fill','on');
            axis equal
            grid on
            axis(xy_lim)
            title(sprintf('Call #%02d',iC));
            colormap(parula(length(cvec)-1));
            caxis([cvec(end), 0])
        end
    end
    
    if plot_opt_all_click
        figure(fig_clicks)
        suptitle(sprintf('%s, all clicks',regexprep(save_fname,'_','\\_')));
        saveSameSize(fig_clicks,'file',...
            fullfile(save_path,sprintf('%s_all_clicks.png',save_fname)),...
            'format','png','renderer','painters');
        close(fig_clicks)
        
        figure(fig_clicks_ofbfix)   % delete out-of-bnd points
        suptitle(sprintf('%s, all clicks',regexprep(save_fname,'_','\\_')));
        saveSameSize(fig_clicks_ofbfix,'file',...
            fullfile(save_path,sprintf('%s_all_clicks_ofbfix.png',save_fname)),...
            'format','png','renderer','painters');
        close(fig_clicks_ofbfix)
    end
    
end
warning on
toc

if save_opt==1
    A.raw_meas = raw_meas;
    A.rotate_data = rotate_data;
    A.shift_tilt_final = shift_tilt_final;
    save(fullfile(save_path,sprintf('%s_all_click_rotated_%s',script_name,bat{iBAT})),'-struct','A');
end

clear A raw_meas rotate_data shift_tilt_final

end