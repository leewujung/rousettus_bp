% 2015 12 03  Merge data from all clicks
% 2015 12 12  Use 35 kHz rotation for other frequencies
% 2015 12 22  Adopt the merging click code for finding multi-freq center of
%             the beampatterns
% 2015 12 24  Move the plotting part to 'compile_info_plot.m'
%             This is only for processing
% 2015 12 24  Separate out the shift and rotation portion and save data
% 2016 04 19  Deal with "out-of-boundary" problem
%             updated shift_rotate_bp.m and bp_fit_ellipse_azel.m
% 2016 07 21  Update the plots to be using map frames
%             Save MAT files for individual clicks
% 2016 10 24  Take out the "out-of-boundary" section and let mic locations
%             "wrap" around back of the bat
%             look for changes by "2016/10/25" in comments


clear
warning off
if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');

    base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/Dropbox/Z_wjlee/projects/rousettus_bp/analysis_results_figs';
else
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

    if strcmp(usrn,'Wu-Jung')
        base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
        save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    else
        base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
    end
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
results_path = 'analysis_results_figs';
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

bat_proc_path = 'proc_output_rousettus_checked';
bat_proc_file_all = dir(fullfile(base_path,bat_proc_path,'rousettus_20150825_*'));

S.base_path = base_path;
S.save_path = save_path;
S.bp_processed_path = bat_proc_path;
S.freq_used_for_rotation = freq_rot;
S.iterative_shift_threshold = it_shift_th;
S.shift_bound_azel = azel_bnd;

% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

S.map.map_projection = map_proj;
S.map.mstruct = mstruct;

% Process all files
for iB = 1:length(bat_proc_file_all)

    tic

    bat_proc_file = bat_proc_file_all(iB).name;
    ss = strsplit(bat_proc_file,'_');
    save_fname = strjoin([script_name,ss(3:4)],'_');

    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(base_path,bat_proc_path,bat_proc_file));
    good_call_idx = find(data.proc.chk_good_call);
    
    S.raw_meas.good_call_idx = good_call_idx;
    
    if plot_opt_all_click
        numrow = ceil(length(good_call_idx)/4);
        fig_clicks = figure;   % shift/tilt
        set(fig_clicks,'Position',[100 100 1400 240*numrow]);
%         fig_clicks_ofbfix = figure;   % shift/tilt and delete out-of-bnd points   % deleted 2016/10/25
%         set(fig_clicks_ofbfix,'Position',[100 100 1400 240*numrow]);
    end
    
    for iC = good_call_idx'

        iC_save = find(iC==good_call_idx);
        fprintf('Call %02d\n',iC);
        
        % Get call info and rotate beampattern using 35 kHz beampattern
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_rot,iC);
        [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
            shift_rotate_bp(az(ch_include_idx),el(ch_include_idx),...
                            call_dB(ch_include_idx),'eckert4',it_shift_th);
        
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
            fig_fit = plot_indiv_click_rotate(Vpass,cvec,mstruct);
            suptitle(sprintf('%s, Call #%02d',regexprep(save_fname,'_','\\_'),iC));
            saveSameSize(fig_fit,'file',...
                         fullfile(save_path,sprintf('%s_c%02d_contour.png',save_fname,iC)),...
                         'format','png','renderer','painters');
            close(fig_fit)
        end

        % For saving MAT files for individual clicks
        S.raw_meas.click_num = iC*ones(length(az),1);
        S.raw_meas.az = az;
        S.raw_meas.el = el;
        S.raw_meas.call_dB = call_dB;
        S.raw_meas.ch_include_idx = ch_include_idx;
        S.raw_meas.click_side = click_side;
        S.raw_meas.click_side_rep = click_side*ones(length(az),1);
        
        S.shift_tilt_final.az = rot_elpctr_tilt.az;
        S.shift_tilt_final.el = rot_elpctr_tilt.el;
%         S.shift_tilt_final.outofbnd_azel_idx = ...      % commented 2016/10/25
%             rot_elpctr_tilt.outofbnd_azel_idx;

        S.raw = raw;
        S.rot_max = rot_max;
        S.rot_elpctr = rot_elpctr;
        S.rot_elpctr_tilt = rot_elpctr_tilt;
        

        % Plot rotated bp using ellipse center
        if plot_opt_all_click
            figure(fig_clicks)
            subplot(numrow,4,find(iC==good_call_idx));
            axesm(mstruct)
%             contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,raw.vq_norm,cvec,'fill','on');
            contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,rot_elpctr_tilt.vq_norm,cvec,'fill','on');  % changed 2016/10/25
            framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
            gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
            tightmap
            axis off
            title(sprintf('Call #%02d',iC));
            colormap(parula(length(cvec)-1));
            caxis([cvec(end), 0])
            
%             figure(fig_clicks_ofbfix)   % delete out-of-bnd points    % deleted 2016/10/25
%             subplot(numrow,4,find(iC==good_call_idx));
%             axesm(mstruct)
%             vq_norm_fix = raw.vq_norm;
%             vq_norm_fix(rot_elpctr_tilt.outofbnd_azel_idx_q) = NaN;
%             contour(rot_elpctr_tilt.xq,rot_elpctr_tilt.yq,...
%                     vq_norm_fix,cvec,'fill','on');
%             framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
%             gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
%             tightmap
%             axis off
%             title(sprintf('Call #%02d',iC));
%             colormap(parula(length(cvec)-1));
%             caxis([cvec(end), 0])
        end

        if save_opt==1
            save(fullfile(save_path,sprintf('%s_c%02d_rotated.mat',save_fname,iC)),...
                 '-struct','S');
        end

    end  % individual click loop
    
    if plot_opt_all_click
        figure(fig_clicks)
        suptitle(sprintf('%s, all clicks',regexprep(save_fname,'_','\\_')));
        saveSameSize(fig_clicks,'file',...
                     fullfile(save_path,sprintf('%s_all_clicks.png',save_fname)),...
                     'format','png','renderer','painters');
        close(fig_clicks)
        
%         figure(fig_clicks_ofbfix)   % delete out-of-bnd points    % deleted 2016/10/25
%         suptitle(sprintf('%s, all clicks',regexprep(save_fname,'_','\\_')));
%         saveSameSize(fig_clicks_ofbfix,'file',...
%                      fullfile(save_path,sprintf('%s_all_clicks_ofbfix.png',save_fname)),...
%                      'format','png','renderer','painters');
%         close(fig_clicks_ofbfix)
    end

    toc
    
end  % each bat_proc file
warning on
%toc

%if save_opt==1
%A.raw_meas = raw_meas;
%A.rotate_data = rotate_data;
%A.shift_tilt_final = shift_tilt_final;
%save(fullfile(save_path,sprintf('%s_all_click_rotated_%s',script_name,bat{iBAT})),'-struct','A');
%end

%clear A raw_meas rotate_data shift_tilt_final

%end