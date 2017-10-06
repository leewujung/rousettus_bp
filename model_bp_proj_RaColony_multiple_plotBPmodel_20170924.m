% 2015 11 24  Use real measurement az/el angles for model beampattern
% 2016 04 19  Use measured mean azimuth width to infer piston model diameter
%             Most processing and plotting codes are from rotate_all_click_20160419.m
% 2016 05 10  Re-process for plots in paper
% 2016 07 21  Add noise at different levels to model
% 2016 07 25  Set flag when shift_rotate_bp error out, see notes for detail
% 2016 08 04  Adapt the code for projecting phased array model using
%             bullethead
% 2016 10 25  Update for version 1025 with out-of-bound points
% 2017 09 21  -- Enable using different beam center point for projection
%             -- Saving info that will allow projection using multiple
%                frequencies in model_composite_RaColony3456_20170921.m
% 2017 09 22  Take out the part to actually carry out simulation
%             and only obtain critical params for making it happen.
%             The actual simulation will happen in model_composite_RaColony3456_20170921.m
%             This modification is done to enable multi-freq simulation
% 2017 09 24  Randomize model bp
% 2017 10 05  Adapt from `model_bp_proj_RaColony_multiple_diffonly_20170924`
%             to plot all BP model used and their beam centers

clear

if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/beampattern_processing');
    
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
    model_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling/';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\beampattern_processing');
    
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
    model_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
end

% Set map projection
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);


% Set params
results_path = 'analysis_results_figs';
rotate_data_path = 'rotate_all_click_2010308';  % rotation on 20170308, note
                                                % the missing '7' in filename
rng_seed = 168;
bpctr_opt = 'ectr';  % 'max' -- max beam energy location
                     % 'top' -- averaged loc for all normalized beam energy>-1
                     % 'ectr' -- center of best-fitting ellipse
it_shift_th = 0.005;
freq_model = 35e3;  % project 35 kHz model and also use 35 kHz for rotation
cgrey = 200*ones(1,3)/255;
cvec = -30:3:0;


% Set save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end


% Set bp model prediction path/files
bp_model_path = 'model_bp_save_multiple_multifreq_20170924';
bp_model_file_all = dir(fullfile(model_base_path,bp_model_path,...
                    sprintf(['model_bp_save_multiple_multifreq_20170924_Ra-colony-rotear-' ...
                    '0.5mm_%dkHz_*.mat'],freq_model/1e3)));
num_bp = length(bp_model_file_all);

for iBP = 1:num_bp
    BP{iBP} = load(fullfile(model_base_path,bp_model_path,bp_model_file_all(iBP).name));

    BP{iBP}.bp_model_path = bp_model_path;
    BP{iBP}.bp_model_file = bp_model_file_all(iBP).name;

    % Load model prediction
    idxsmall = BP{iBP}.pp_plot<-30;
    BP{iBP}.pp_plot(idxsmall) = -30;
    idxnotnan = ~isnan(BP{iBP}.pp_plot);
    [~,BP{iBP}.vq_norm,BP{iBP}.azq,BP{iBP}.elq] = ...
        interp_bp(BP{iBP}.az(idxnotnan)/180*pi,...
                  BP{iBP}.el(idxnotnan)/180*pi,BP{iBP}.pp_plot(idxnotnan),'rbf');
    BP{iBP}.azq = BP{iBP}.azq/pi*180;
    BP{iBP}.elq = BP{iBP}.elq/pi*180;

    % Find beam center for model
    % --- max beam energy location
    xx = BP{iBP}.vq_norm(:);
    [~,mm_idx] = max(xx);
    BP{iBP}.max.el = BP{iBP}.elq(mm_idx);
    BP{iBP}.max.az = BP{iBP}.azq(mm_idx);
    % --- averaged location of normalized beam energy >-1
    xx(isnan(xx)) = -inf;
    [~,sort_idx] = sort(xx,'descend');
    ii = xx(sort_idx)>-1;
    BP{iBP}.top.el = mean(BP{iBP}.elq(sort_idx(ii)));
    BP{iBP}.top.az = mean(BP{iBP}.azq(sort_idx(ii)));
    % --- fit ellipse
    % ------- BP{iBP}.el = bem_results.theta
    % ------- BP{iBP}.az = bem_results.phi
    [raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(BP{iBP}.az,BP{iBP}.el,BP{iBP}.pp,map_proj,0.005);
    [el_ectr,az_ectr] = minvtran(mstruct,rot_max.E.x0,rot_max.E.y0);  % inverse map projection
    [BP{iBP}.ectr.el,BP{iBP}.ectr.az] = rotatem(el_ectr,az_ectr,...
                                                [BP{iBP}.max.el,BP{iBP}.max.az],...
                                                'inverse','degrees');

    % Plot BP model
    figure
    plot_bp_simple(subplot(1,1,1),...
                   -BP{iBP}.az,BP{iBP}.el,BP{iBP}.pp_plot,map_proj);
    hold on
    plotm(BP{iBP}.max.el,-BP{iBP}.max.az,'rx','markersize',8,'linewidth',2);
    plotm(BP{iBP}.top.el,-BP{iBP}.top.az,'r^','markersize',8,'linewidth',2);
    plotm(BP{iBP}.ectr.el,-BP{iBP}.ectr.az,'ro','markersize',8,'linewidth',2);
    title('model: left');
    gridm('gcolor',cgrey,'glinestyle','-');
    framem('fedgecolor',cgrey);
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    save_fname = sprintf('%s_%s',script_name,bp_model_file_all(iBP).name);
    saveSameSize_res(gcf,120,'file',fullfile(save_path,[save_fname,'.png']),...
                 'format','png','renderer','painters');
    epswrite(fullfile(save_path,[save_fname,'.eps']));
    
end


