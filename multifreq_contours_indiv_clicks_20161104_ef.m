% 2015 12 27  Plot contours of the same frequency across all clicks together

% clear
%
% usrn = getenv('username');
% addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
% addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
%
% plot_opt = 1;
% save_opt = 1;
%
% % Load compiled rotated data
% base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\analysis_results_figs';
% compile_path = '20151224_all_click_rotate';
% compile_file = 'all_click_rotated_39184.mat';
% load(fullfile(base_path,compile_path,compile_file));

% Fix path with current username
% if ~strcmp(usrn,'wlee76')
%     if strfind(base_path, 'wlee76')
%         base_path = regexprep(base_path,'wlee76',usrn);
%     end
% end

% A.base_path = base_path;
% A.compile_path = compile_path;
% A.compile_file = compile_file;

% % Path/filename for saving data and figures
% save_header = 'eptesicus_20150824_LB53';
% save_path = fullfile(base_path,'20151228_indiv_click_contour');  % set path for saving files
% if ~exist(save_path,'dir')
%     mkdir(save_path);
% end
%
% Frequency to be plotted
freq_wanted = (30:5:60)*1e3;
num_freq = length(freq_wanted);

% A.freq_wanted  = freq_wanted;

% processed_files

% warning off
% Get contours from all freq
% c3db_xy_all = cell(1,1);
contour_sm_len = 10;
colorset = jet(num_freq);


% for iB=1:length(processed_files)
bat_proc_file = 'eptesicus_20150824_LB53_18_mic_data_bp_proc.mat';
data_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\proc_output_eptesicus_new';
ss = strsplit(bat_proc_file,'_');
save_fname = strjoin(ss(3:4),'-');
fprintf('File: %s\n',bat_proc_file);

fig_mf_indiv = figure;   % plot individual multifreq contour
axesm eckert4
axis off
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');

data = load(fullfile(data_path,bat_proc_file));
good_call_idx = find(data.proc.chk_good_call);
iB=1;
iC = 20%good_call_idx'
iC_save = find(iC==good_call_idx);
fprintf('Call %02d\n',iC);

for iF=1:num_freq
    % Get data and good mic index
    [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted(iF),iC);
    call_dB_norm = call_dB - max(call_dB);  % normalize
    [~,vq_norm,azq,elq] = interp_bp(az,el,call_dB,'rbf');  % use the first frequency data for finding ellipse center
    azq=azq/pi*180;
    elq=elq/pi*180;
    % Get -3dB contour
    %             M = rotate_data{iB}(iC_save).rot_elpctr_tilt;
    idx_notnan = ~isnan(azq(:));  % index of non-NaN data
    %             azq = M.azq(idx_notnan);
    %             elq = M.elq(idx_notnan);
    azq = min(azq(:)):max(azq(:));
    elq = min(elq(:)):max(elq(:));
    [azq,elq] = meshgrid(azq,elq);
    %             [azq,elq,vq] = griddata(M.azq(idx_notnan),M.elq(idx_notnan),vq_norm(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq
    [azq,elq,vq] = griddata(azq(idx_notnan),elq(idx_notnan),vq_norm(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq
    
    [~,c_main_nan] = get_main_contour(vq,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
    [c3db_xy(:,1),c3db_xy(:,2)] = mfwdtran(mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]
    %             c3db_xy_all{iB}{iF,iC_save} = c3db_xy;
    
    
    figure(fig_mf_indiv)
    xy_sm(:,1) = smooth(c3db_xy(:,1),contour_sm_len);
    xy_sm(:,2) = smooth(c3db_xy(:,2),contour_sm_len);
    xy_sm(isnan(c3db_xy(:,1)),:) = NaN;
    plot(xy_sm(:,1),xy_sm(:,2),'linewidth',2,'color',colorset(iF,:));
    hold on
    clear xy_sm
    clear c3db_xy
    
end

colormap(jet(num_freq))
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
    'TickLabels',{num2str(freq_wanted'/1e3)},'location','southoutside');
title(sprintf('%s, Call #%02d',regexprep(save_fname,'_','\\_'),iC));

