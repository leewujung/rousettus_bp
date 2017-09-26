% 2017 09 25  Compare the beam center location for data and model, use all
%             three possible defitions of beam center


clear

usrn = getenv('username');
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    data_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    data_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp';
end

results_path = 'analysis_results_figs';
model_path = 'model_multifreq_indiv_bpctr_cntr_20170924_std1.0';
model_file = 'model_multifreq_indiv_bpctr_cntr_20170924_std1.0_results.mat';
data_path = 'multifreq_indiv_bpctr_cntr_20170924';
data_file = 'multifreq_indiv_bpctr_cntr_20170924_batall_results.mat';

model = load(fullfile(data_base_path,results_path,model_path,model_file));
data = load(fullfile(data_base_path,results_path,data_path,data_file));

[~,model_freqI] = sort(model.freq.all);
freq_idx_wanted = 2:length(model.freq.all)-1;


% Set save folder
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(data_base_path,results_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

save_fname = script_name;


% Plot mean beam center az locations
fig_bpctr_stat = figure('position',[180,100,880,850]);

% ==================================================
% Left clicks, bpctr_opt = max
subplot(321)
plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_left],'o-')
hold on
plot([data.bpctr_az.mean(freq_idx_wanted).max_left],'o-')
legend('Model','Data','location','best');
title('Left, beam center = max')
grid; ylim([-15,15]);
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% Right clicks, bpctr_opt = max
subplot(322)
plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_right],'o-')
hold on
plot([data.bpctr_az.mean(freq_idx_wanted).max_right],'o-')
legend('Model','Data','location','best')
title('Right, beam center = max')
grid; ylim([-15,15]);
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% ==================================================
% Left clicks, bpctr_opt = top
subplot(323)
plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_left],'o-')
hold on
plot([data.bpctr_az.mean(freq_idx_wanted).top_left],'o-')
legend('Model','Data','location','best')
title('Left, beam center = top')
grid; ylim([-15,15]);
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% Right clicks, bpctr_opt = top
subplot(324)
plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_right],'o-')
hold on
plot([data.bpctr_az.mean(freq_idx_wanted).top_right],'o-')
legend('Model','Data','location','best')
title('Right, beam center = top')
grid; ylim([-15,15]);
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% ==================================================
% Left clicks, bpctr_opt = ectr
subplot(325)
plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_left],'o-')
hold on
plot([data.bpctr_az.mean(freq_idx_wanted).ectr_left],'o-')
legend('Model','Data','location','best')
title('Left, beam center = ectr')
grid; ylim([-15,15]);
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% Right clicks, bpctr_opt = ectr
subplot(326)
plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_right],'o-')
hold on
plot([data.bpctr_az.mean(freq_idx_wanted).ectr_right],'o-')
legend('Model','Data','location','best')
title('Right, beam center = ectr')
grid; ylim([-15,15]);
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')




saveas(fig_bpctr_stat,fullfile(save_path,sprintf('%s.fig',save_fname)),'fig');
saveSameSize_res(fig_bpctr_stat,150,'file',...
                 fullfile(save_path,sprintf('%s.png',save_fname)),...
                 'format','png','renderer','painters');
epswrite(fullfile(save_path,sprintf('%s.eps',save_fname)));
