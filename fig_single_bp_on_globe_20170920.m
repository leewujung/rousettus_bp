% 2016 04 15  Plot single beam pattern on globe, eckert4 projection
%             Code adopted from plot_bp_on_globe.m
% 2016 05 07  Plot for NIFTI poster
% 2016 05 07  Plot for ASA talk
% 2016 07 21  Plot for paper
% 2017 09 20  -- Add center of the beam on the figure
%             -- Check location of the center of best-fitting ellipse

clear
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');

    data_base_path = '~/internal_2tb/Dropbox/0_ANALYSIS/bp_processing';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\export_fig-master');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    
    data_base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end

% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));

results_path = 'analysis_results_figs';
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

bat_proc_path = 'proc_output_rousettus_new';
bat_proc_file = 'rousettus_20150825_36134_02_mic_data_bp_proc.mat';

%freq_wanted = 35*1e3;
freq_wanted = [20:5:55]*1e3;
map_proj_opt = 'eckert4';
% map_proj_opt = 'ortho';
interp_opt = 'rb_rbf';

call_num = 19:20;

data = load(fullfile(data_base_path,bat_proc_path,bat_proc_file));
ss = strsplit(bat_proc_file,'_');

for iF=1:length(freq_wanted)
    for iC=call_num
        ss_title = sprintf('bat %s, trial %s, click #%02d, freq %d kHz',ss{3},ss{4},iC,freq_wanted(iF)/1e3);

        % Load data from a particular freq
        [call_dB,az,el,ch_good_idx] = get_call_azel_dB_data(data,freq_wanted(iF),iC);

        % Take into account only good mics
        az = az(ch_good_idx);
        el = el(ch_good_idx);
        call_dB = call_dB(ch_good_idx);
                
        % Normalize and interpolation
        maxref = max(call_dB);
        call_dB_norm = call_dB-maxref;
        [azq,elq] = meshgrid(min(az):pi/180:max(az),min(el):pi/180:max(el));
        if strcmp(interp_opt,'rb_natural')  % natural neighbor interpolation
            vq = griddata(az,el,call_dB(angle_notnanidx),azq,elq,'natural');
        elseif strcmp(interp_opt,'rb_rbf')  % radial basis function interpolation
            vq = rbfinterp([azq(:)';elq(:)'],rbfcreate([az(:)';el(:)'],call_dB(:)','RBFFunction','multiquadrics'));
            vq = reshape(vq,size(azq));
        end
        vq_norm = vq-maxref;
        
        % Find indices within measured polygon
        k = boundary(az,el,0);  % outer boundary of all measured points
        [in,on] = inpolygon(azq,elq,az(k),el(k));
        in_smpl_poly = in|on;
        clear in on
        vq(~in_smpl_poly) = NaN;
        vq_norm(~in_smpl_poly) = NaN;

        % Plot interpolated beampattern with eckert4 projection
        elqm = elq;
        elqm(isnan(vq_norm)) = NaN;
        azqm = azq;
        azqm(isnan(vq_norm)) = NaN;

        % Find beam center (=max beam energy point)
        xx = vq_norm(:);
        [~,vq_norm_max_idx] = max(xx);
        vq_norm_max_el_loc = elq(vq_norm_max_idx);
        vq_norm_max_az_loc = azq(vq_norm_max_idx);
        xx(isnan(xx)) = -inf;
        [~,sort_idx] = sort(xx,'descend');
        ii = xx(sort_idx)>-1;
        vq_norm_top_el_loc = mean(elq(sort_idx(ii)));
        vq_norm_top_az_loc = mean(azq(sort_idx(ii)));
        
        % Fit ellipse
        map_proj = 'eckert4';
        mstruct = defaultm(map_proj);
        mstruct = defaultm(mstruct);
        [click_side,raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
            shift_rotate_bp(az(ch_good_idx),el(ch_good_idx),...
                            call_dB(ch_good_idx),map_proj,0.005);

        % Project and rotate the ellipse center back to the map/beam coordinate
        [el_ectr,az_ectr] = minvtran(mstruct,rot_max.E.x0,rot_max.E.y0);  % inverse map projection
        [el_ectr_r,az_ectr_r] = rotatem(el_ectr,az_ectr,...
                                        [vq_norm_max_el_loc,vq_norm_max_az_loc]/pi*180,...
                                        'inverse','degrees');
        
        %     vq_norm_min = min(vq_norm(:));  % routine used in bp GUI
        %     contour_vec = -3:-3:(floor(vq_norm_min/3)-1)*3;
        vq_norm_min = -27;
        contour_vec = 0:-3:(floor(vq_norm_min/3)-1)*3;
        cvec_min_idx = find(contour_vec-vq_norm_min<0,1,'first');

        % Plot bp
        fig_bp = figure;
        axesm(map_proj_opt);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        framem('fedgecolor',190*ones(1,3)/255,'flonlimit',[-180 180]);
        axis off
        contourfm(elq/pi*180,azq/pi*180,vq_norm,contour_vec(2:cvec_min_idx),...
            'fill','on','linecolor','w');  % don't plot 0 dB contour
        hold on
        
        % location of mics
        plotm(el/pi*180,az/pi*180,'kx','markersize',8,'linewidth',2);
        
        % location of max beam energy
        plotm(vq_norm_max_el_loc/pi*180,vq_norm_max_az_loc/pi*180,...
              'rx','markersize',8,'linewidth',2)

        % location of max beam energy
        plotm(vq_norm_top_el_loc/pi*180,vq_norm_top_az_loc/pi*180,...
              'r^','markersize',8,'linewidth',2)

        % location of center of best-fitting ellipse
        plotm(el_ectr_r,az_ectr_r,'ro','markersize',8,'linewidth',2)
        
%         textm(el/pi*180,az/pi*180,num2str(mic_num(angle_notnanidx)'),...
%             'horizontalalignment','center','fontsize',14);
        colorbar('southoutside','ticks',fliplr(contour_vec(1:cvec_min_idx)));
        colormap(parula(cvec_min_idx-1));
        caxis([contour_vec(cvec_min_idx) 0]);
        tightmap
        title(ss_title)
        set(fig_bp,'color','w');
        
        ss_file = sprintf('%s_%s_%s_c%02d_f%dkHz_%s_%s_mic',...
            script_name,ss{3},ss{4},iC,freq_wanted(iF)/1e3,map_proj_opt,interp_opt(4:end));

        saveas(fig_bp,fullfile(save_path,[ss_file,'.fig']),'fig');
        saveSameSize_300(fig_bp,'file',fullfile(save_path,ss_file),...
            'format','png','renderer','painters');
        epswrite(fullfile(save_path,[ss_file,'.eps']));
       
    end
end

