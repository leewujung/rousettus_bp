% 2016 08 08  Assemble composite clicks from simulated mic receptions
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 09 21  Allow assembling composite clicks from multiple frequencies

clear

usrn = getenv('username');
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end

save_opt = 1;

results_path = 'analysis_results_figs';
simu_data_path = 'model_bp_proj_RaColony_20170921_std1.0';
ss = strsplit(simu_data_path,'_');
model_shape = ss{end-2};
model_calc_date = ss{end-1};
noise_std = ss{end};
simu_file = dir(fullfile(data_base_path,results_path,simu_data_path,'*.mat'));

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

A.param.data_base_path = data_base_path;
A.param.simu_data_path = simu_data_path;
A.param.simu_data_file = simu_file;

% Composite click params
threshold = 0;  % only use averaged results if number of measurement > threshold
binsize = 10;  % bin size in azimuth and elevation

A.param.composite_threshold = threshold;
A.param.composite_binsize = binsize;

% Load 1 file to get mtx sizes
D = load(fullfile(data_base_path,results_path,simu_data_path,simu_file(1).name));  % rotated data
num_ch = length(D.model_rot.raw.az);

A.param.num_ch = num_ch;
A.param.map = D.map;

freq_model_all = [25:10:55]*1e3;
noise_std = 1;

for iF=1:length(freq_model_all)

    freq_model = freq_model_all(iF);

    % Set save_fname
    ss_save = strsplit(D.BP.bp_model_file,'_');
    tongue_loc_str = strtok(strjoin({ss_save{end-2:end}},'_'),'.');

    save_fname = sprintf('%s_noise%2.1f_%dkHz_%s_bin%d_th%d',...
                         script_name,noise_std,freq_model/1e3,...
                         tongue_loc_str,binsize,threshold);


    click_side_all = nan(length(simu_file),num_ch);
    az_elpctr_tilt_all = nan(length(simu_file),num_ch);
    el_elpctr_tilt_all = nan(length(simu_file),num_ch);
    call_dB_all = nan(length(simu_file),num_ch);
    for iS=1:length(simu_file)

        D = load(fullfile(data_base_path,results_path,simu_data_path, ...
                          simu_file(iS).name));  % load simulation results
        for iFO = 1:length(freq_other)
            bp_model_file_other = ...
                sprintf(['model_bp_save_20161009_multifreq_Ra-colony-rotear-' ...
                         '0.5mm_%dkHz_x029_y000_z-05.mat'],freq_other(iFO)/1e3);
            BP_other = load(fullfile(model_base_path,bp_model_path,bp_model_file_other));
            idxnotnan_other = ~isnan(BP_other.pp_plot);
            vq_mic_other(:,iFO) = rbfinterp([az(:)';el(:)'],...
                                            rbfcreate([az_model_rot(idxnotnan_other)';...
                                el_model_rot(idxnotnan_other)'],...
                                                      BP_other.pp_plot(idxnotnan_other)',...
                                                      'RBFFunction','multiquadrics'));
        end
        vq_mic_other = vq_mic_other+randn(size(vq_mic_other))*noise_std+noise_mean;


        % Determine right/left click
        [~,mmidx]=max(D.raw.vq_norm(:));
        if D.raw.azq(mmidx)<0
            click_side = 0;  % left: click_side = 0
        else
            click_side = 1;  % right: click_side = 1
        end

        % Merge all data
        click_side_all(iS,:) = click_side*ones(1,length(num_ch));
        if isempty(D.rot_elpctr_tilt)
            az_elpctr_tilt_all(iS,:) = nan(size(D.raw.az));
            el_elpctr_tilt_all(iS,:) = nan(size(D.raw.az));
            call_dB_all(iS,:) = nan(size(D.raw.az));
        else
            az_elpctr_tilt_all(iS,:) = D.rot_elpctr_tilt.az;
            el_elpctr_tilt_all(iS,:) = D.rot_elpctr_tilt.el;
            call_dB_all(iS,:) = D.rot_elpctr_tilt.call_dB;
        end

    end
    A.click_side_all = click_side_all;
    A.az_elpctr_tilt_all = az_elpctr_tilt_all;
    A.el_elpctr_tilt_all = el_elpctr_tilt_all;
    A.call_dB_all = call_dB_all;




    % All data points for a single freq
    call_dB_merge = reshape(call_dB_all,[],1);
    click_side_merge = reshape(click_side_all,[],1);
    az_elpctr_tilt_merge = reshape(az_elpctr_tilt_all,[],1);
    el_elpctr_tilt_merge = reshape(el_elpctr_tilt_all,[],1);

    idx_right_good = find(click_side_merge==1 & ~isnan(call_dB_merge));
    call_dB_right = call_dB_merge(idx_right_good);
    az_right = az_elpctr_tilt_all(idx_right_good);
    el_right = el_elpctr_tilt_all(idx_right_good);

    idx_left_good = find(click_side_merge==0 & ~isnan(call_dB_merge));
    call_dB_left = call_dB_merge(idx_left_good);
    az_left = az_elpctr_tilt_all(idx_left_good);
    el_left = el_elpctr_tilt_all(idx_left_good);

    % Interpolation using averaged data
    [interp_avg_right,bin_avg_right] = average_call(az_right,el_right,call_dB_right,...
                                                    binsize,'eckert4',threshold);
    [interp_avg_left,bin_avg_left] = average_call(az_left,el_left,call_dB_left,...
                                                  binsize,'eckert4',threshold);

    % Load model bp
    BP = load(fullfile(D.bp_model_path,D.bp_model_file));
    idxnotnan = ~isnan(BP.pp_plot);
    [~,BP.vq_norm,BP.azq,BP.elq] = ...
        interp_bp(BP.az(idxnotnan)/180*pi,BP.el(idxnotnan)/180*pi,BP.pp_plot(idxnotnan),'rbf');
    BP.azq = BP.azq/pi*180;
    BP.elq = BP.elq/pi*180;

    % Plot params
    suptitle_text = sprintf('%s, %s, %s',model_shape,noise_std,...
                            regexprep(D.bp_model_file,'_','\\_'));
    cgrey = 200*ones(1,3)/255;
    cvec=-30:3:0;

    % Plot all scatter samples
    fig_scatter = figure('position',[200,200,1200,450]);

    subplot(121)  % left click
    axesm(D.map.map_projection);
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    scatterm(el_left,az_left,15,call_dB_left,'filled')
    title('Left click');
    tightmap
    axis off
    caxis(cvec([1 end]))
    colorbar('location','southoutside');

    subplot(122)  % right click
    axesm(D.map.map_projection);
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    scatterm(el_right,az_right,15,call_dB_right,'filled')
    title('Right click');
    tightmap
    axis off
    caxis(cvec([1 end]))
    colorbar('location','southoutside');


    % Plot model bp and reconstructed bp for comparison
    mstruct = defaultm(D.map.map_projection);
    mstruct = defaultm(mstruct);
    [BP.xq,BP.yq] = mfwdtran(mstruct,BP.elq,BP.azq);
    BP.vq_norm(BP.vq_norm<min(cvec)) = min(cvec);


    fig_bp_cmp = figure('position',[200,60,1200,900]);

    subplot(221)    % model bp: left
    axesm(D.map.map_projection);
    contourf(-BP.xq,BP.yq,BP.vq_norm,cvec,'w');
    % contourm(BP.elq,-BP.azq,BP.vq_norm,cvec,'fill','on');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    title('Model bp: left');
    tightmap
    axis off
    colormap(parula(length(cvec)-1))
    caxis(cvec([1 end]))
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    subplot(222)    % model bp: right
    axesm(D.map.map_projection);
    contourf(BP.xq,BP.yq,BP.vq_norm,cvec,'w');
    % contourm(BP.elq,BP.azq,BP.vq_norm,cvec,'fill','on');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    title('Model bp: right');
    tightmap
    axis off
    colormap(parula(length(cvec)-1))
    caxis(cvec([1 end]))
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    subplot(223)   % reconstructed: left
    axesm(D.map.map_projection);
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    contourfm(interp_avg_left.elq_avg,interp_avg_left.azq_avg,interp_avg_left.vq_norm_avg,cvec,'w');
    title('Reconstructed bp: left');
    tightmap
    axis off
    colormap(parula(length(cvec)-1))
    caxis(cvec([1 end]))
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    subplot(224)   % reconstructed: right
    axesm(D.map.map_projection);
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    contourfm(interp_avg_right.elq_avg,interp_avg_right.azq_avg,interp_avg_right.vq_norm_avg,cvec,'w');
    title('Reconstructed bp: right');
    tightmap
    axis off
    colormap(parula(length(cvec)-1))
    caxis(cvec([1 end]))
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    mtit(suptitle_text)

    % Save figure and mat file

    figure(fig_scatter)
    saveSameSize(fig_scatter,'file',fullfile(save_path,[save_fname,'_scatter.png']),...
                 'format','png','renderer','painters');
    epswrite(fullfile(save_path,[save_fname,'_scatter.eps']));

    figure(fig_bp_cmp)
    saveSameSize(fig_bp_cmp,'file',fullfile(save_path,[save_fname,'_avg_bp.png']),...
                 'format','png','renderer','painters');
    epswrite(fullfile(save_path,[save_fname,'_avg_bp.eps']));

    save(fullfile(save_path,[save_fname,'.mat']),'-struct','A');

end % loop through all freq