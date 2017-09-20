% 2015 12 27  Plot contours of the same frequency across all clicks together
% 2016 04 19  Revisit code, update data format from rotate_all_click.m
% 2016 05 07  Plot for NIFTI poster

clear

usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

save_opt = 1;
plot_contour_indiv = 1;
plot_contour_all = 0;

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_root = 'analysis_results_figs';
save_path = fullfile(base_path,save_root,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

data_path = 'rotate_all_click_20160419';


% Params
freq_wanted = (20:5:50)*1e3;
num_freq = length(freq_wanted);
contour_sm_len = 10;
colorset = jet(num_freq);

bat = {'3bat','36134','34271','39184'};

for iBAT=2%:length(bat)
    
    data_file = ['rotate_all_click_20160419_all_click_rotated_',bat{iBAT},'.mat'];
    save_fname = ['indiv_click_contour_',bat{iBAT},'.mat'];
    
    A.data_file = data_file;
    A.base_path = base_path;
    A.data_path = data_path;
    A.freq_wanted  = freq_wanted;
    A.contour_smooth_len_for_plot = contour_sm_len;  % for plotting only, not in stored data
    
    D = load(fullfile(base_path,save_root,data_path,data_file));
    
    warning off
    
    % Get contours from all freq
    c3db_xy_all = cell(length(D.rotate_data),1);
    for iB=1%:length(D.rotate_data)
        ss = strsplit(D.bp_processed_file_all(iB).name,'_');   % new format
        save_fname = strjoin([script_name,ss(3:4)],'_');
        fprintf('File: %s\n',D.bp_processed_file_all(iB).name);
        
        data = load(fullfile(base_path,D.bp_processed_path,D.bp_processed_file_all(iB).name));
        good_call_idx = find(data.proc.chk_good_call);
        
        for iC = 19:20 %good_call_idx'
            iC_save = find(iC==good_call_idx);
            fprintf('Call %02d\n',iC);
            
            if plot_contour_indiv
                fig_mf_indiv = figure;   % plot individual multifreq contour
                axesm eckert4
                axis off
                framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
                gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
            end
            
            for iF=1:num_freq
                % Get data and good mic index
                [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted(iF),iC);
                call_dB_norm = call_dB - max(call_dB);  % normalize
                [~,vq_norm,~,~] = interp_bp(az,el,call_dB,'rbf');  % use the first frequency data for finding ellipse center
                
                % Get -3dB contour: rotated, for plotting all contours together
                M = D.rotate_data{iB}(iC_save).rot_elpctr_tilt;
                idx_notnan = ~isnan(M.azq) & ~M.outofbnd_azel_idx_q & ~isnan(M.vq_norm);  % index of non-NaN data
                azq = M.azq(idx_notnan);
                elq = M.elq(idx_notnan);
                azq = min(azq(:)):max(azq(:));
                elq = min(elq(:)):max(elq(:));
                [azq,elq] = meshgrid(azq,elq);
                [azq,elq,vq] = griddata(M.azq(idx_notnan),M.elq(idx_notnan),vq_norm(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq
                
                [~,c_main_nan] = get_main_contour(vq,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
                [c3db_xy(:,1),c3db_xy(:,2)] = mfwdtran(D.map.mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]
                c3db_xy_all{iB}{iF,iC_save} = c3db_xy;
                clear c3db_xy
                
                % Plot
                if plot_contour_indiv
                    % Get -3dB contour: raw, for plotting individual clicks
                    M = D.rotate_data{iB}(iC_save).raw;
                    idx_notnan = ~isnan(M.azq) & ~isnan(M.vq_norm);  % index of non-NaN data
                    azq = M.azq(idx_notnan);
                    elq = M.elq(idx_notnan);
                    azq = min(azq(:)):max(azq(:));
                    elq = min(elq(:)):max(elq(:));
                    [azq,elq] = meshgrid(azq,elq);
                    [azq,elq,vq] = griddata(M.azq(idx_notnan),M.elq(idx_notnan),vq_norm(idx_notnan),azq,elq);  % griddata to create regular orthogonal azq/elq
                    
                    [~,c_main_nan] = get_main_contour(vq,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
                    [c3db_xy(:,1),c3db_xy(:,2)] = mfwdtran(D.map.mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]
                    
                    figure(fig_mf_indiv)
                    xy_sm(:,1) = smooth(c3db_xy(:,1),contour_sm_len);
                    xy_sm(:,2) = smooth(c3db_xy(:,2),contour_sm_len);
                    xy_sm(isnan(c3db_xy(:,1)),:) = NaN;
                    plot(xy_sm(:,1),xy_sm(:,2),'linewidth',3,'color',colorset(iF,:));
                    hold on
                    clear xy_sm
                end
                
                clear c3db_xy
            end
            
            if plot_contour_indiv
                figure(fig_mf_indiv)
                tightmap
                colormap(jet(num_freq))
                colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
                    'TickLabels',{num2str(freq_wanted'/1e3)},'location','southoutside');
                title(sprintf('%s, Call #%02d',regexprep(save_fname,'_','\\_'),iC));
                saveSameSize(fig_mf_indiv,'file',...
                    fullfile(save_path,sprintf('%s_c%02d.png',save_fname,iC)),...
                    'format','png','renderer','painters');
                saveas(fig_mf_indiv,...
                    fullfile(save_path,sprintf('%s_c%02d.fig',save_fname,iC)),'fig');
                close(fig_mf_indiv)
            end
        end
    end
    A.contour_3db_xy = c3db_xy_all;
    
    save_fname_all_click = strjoin([script_name,ss(3)],'_');
    if save_opt==1
        save(fullfile(save_path,sprintf('%s_all click.mat',save_fname_all_click)),'-struct','A');
    end

end % loop through each bat
