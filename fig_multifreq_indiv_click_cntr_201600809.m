% 2015 12 27  Plot contours of the same frequency across all clicks together
% 2016 04 19  Revisit code, update data format from rotate_all_click.m
% 2016 05 07  Plot for NIFTI poster
% 2016 08 09  Plot for paper;
%             revised from multifreq_indiv_click_cntr_extract_20160507.m

clear

usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    data_base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
else
    data_base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\proc_output_rousettus_checked'];
    save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
results_path = 'analysis_results_figs';
save_path = fullfile(save_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

data_file = 'rousettus_20150825_36134_02_mic_data_bp_proc.mat';
ss = strsplit(data_file,'_');
bat_num = ss{3};
trial_num = str2double(ss{4});

% Params
freq_wanted = (25:5:55)*1e3;
num_freq = length(freq_wanted);
contour_sm_len = 10;
colorset = jet(num_freq);

% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

warning off

% Load data
data = load(fullfile(data_base_path,data_file));
good_call_idx = find(data.proc.chk_good_call);

cgrey = 200*ones(1,3)/255;
for iC = 19:20 %good_call_idx'
    iC_save = find(iC==good_call_idx);
    fprintf('Call %02d\n',iC);
    
    fig_mf_indiv = figure;   % plot individual multifreq contour
    axesm eckert4
    axis off
    framem('fedgecolor',cgrey,'flonlimit',[-180 180]);
    gridm('gcolor',cgrey,'glinestyle','-');
    
    for iF=1:num_freq
        % Get data and good mic index
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted(iF),iC);
        call_dB_norm = call_dB - max(call_dB);  % normalize
        [~,vq_norm,azq,elq] = interp_bp(az,el,call_dB,'rbf');  % use the first frequency data for finding ellipse center
        azq = azq/pi*180;
        elq = elq/pi*180;
        idx_notnan = ~isnan(azq) & ~isnan(vq_norm);  % index of non-NaN data

        [azq,elq,vq] = griddata(azq(idx_notnan),elq(idx_notnan),vq_norm(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq
        
        [~,c_main_nan] = get_main_contour(vq,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
        [c3db_xy(:,1),c3db_xy(:,2)] = mfwdtran(mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]
        
        figure(fig_mf_indiv)
        xy_sm(:,1) = smooth(c3db_xy(:,1),contour_sm_len);
        xy_sm(:,2) = smooth(c3db_xy(:,2),contour_sm_len);
        xy_sm(isnan(c3db_xy(:,1)),:) = NaN;
        plot(xy_sm(:,1),xy_sm(:,2),'linewidth',3,'color',colorset(iF,:));
        hold on

        clear xy_sm
        clear c3db_xy
    end
    
    figure(fig_mf_indiv)
    tightmap
    colormap(jet(num_freq))
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(freq_wanted'/1e3)},'location','southoutside');
    title(sprintf('Bat %s, Trial %02d, Call #%02d',bat_num,trial_num,iC));
    
    save_fname = sprintf('%s_%s_%02d_c%02d_mf_cntr',script_name,bat_num,trial_num,iC);
    saveSameSize(fig_mf_indiv,'file',...
        fullfile(save_path,sprintf('%s_c%02d.png',save_fname,iC)),...
        'format','png','renderer','painters');
    saveas(fig_mf_indiv,...
        fullfile(save_path,sprintf('%s_c%02d.fig',save_fname,iC)),'fig');
    close(fig_mf_indiv)
end


