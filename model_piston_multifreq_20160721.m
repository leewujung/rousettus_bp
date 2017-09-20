% 2016 07 21  Model 2D beampattern of circular piston

clear

if isunix
    addpath('~/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/Dropbox/0_CODE/beampattern_processing');

    base_path = '~/Dropbox/Z_wjlee/projects/rousettus_bp/analysis_results_figs';
else
    usrn = getenv('username');
    if strcmp(usrn,'Wu-Jung')
        addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
        addpath('F:\Dropbox\0_CODE\MATLAB\EllipseDirectFit');
        addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
        addpath('F:\Dropbox\0_CODE\beampattern_processing');
    else
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\EllipseDirectFit']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
        addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_processing']);
    end
    
    % Set base path
    if strcmp(usrn,'Wu-Jung')
        base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
    else
        base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
    end

end

noise_mean = 0;  % added noise profile [dB]
noise_std = 1;   % [dB]

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_root = 'analysis_results_figs';
save_path = fullfile(base_path,save_root,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end


% Set bp param
a = 5e-3;  % [m] aperture radius
c = 344;  % sound speed [m/s]
freq_all = [25:5:55]*1e3;
k_all = 2*pi*freq_all/c;  % wavenumber
colorset = jet(length(k_all));

% Calc bp and find -3dB points
pol_angle = (0:0.1:90)/180*pi;
mic_dB = 20*log10(abs(2*besselj(1,k_all'*a*sin(pol_angle))./...
    (k_all'*a*sin(pol_angle))));
[~,bp_idx] = min(abs(mic_dB-(-3)),[],2);

% Plot -3 dB contours on map coord
fig_mf = figure;
axesm('eckert4');
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-90 90]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
axis off
for iK=1:length(k_all)
    x = pol_angle(bp_idx(iK))/pi*180*cos(0:0.01*pi:2*pi);
    y = pol_angle(bp_idx(iK))/pi*180*sin(0:0.01*pi:2*pi);
    
    figure(fig_mf)
    plotm(x,y,'linewidth',2,'color',colorset(iK,:));
    hold on
end
hold off
tightmap
title(sprintf('Circular piston, radius %2.2f mm, eckert4',a*1e3),'fontsize',14)

save_fname = sprintf('circ_piston_r%2.2fmm_eckert4',a*1e3);

saveas(fig_mf,...
    fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig_mf,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');


fig_mf = figure;
axesm('ortho');
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-90 90]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
axis off
for iK=1:length(k_all)
    x = pol_angle(bp_idx(iK))/pi*180*cos(0:0.01*pi:2*pi);
    y = pol_angle(bp_idx(iK))/pi*180*sin(0:0.01*pi:2*pi);
    
    figure(fig_mf)
    plotm(x,y,'linewidth',2,'color',colorset(iK,:));
    hold on
end
hold off
tightmap
title(sprintf('Circular piston, radius %2.2f mm, ortho',a*1e3),'fontsize',14)

save_fname = sprintf('circ_piston_r%2.2fmm_ortho',a*1e3);

saveas(fig_mf,...
    fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig_mf,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');