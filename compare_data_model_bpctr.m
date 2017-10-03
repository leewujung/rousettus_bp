% 2017 09 25  Compare the beam center location for data and model, use all
%             three possible defitions of beam center
% 2017 09 26  MANOVA to see if model and data bpctr can be consider "the
%             same" as a whole across multiple frequencies

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

n = 300;
plot_type = 0; % 0-'meanonly', 1-'withstd'

%% Plot mean beam center az locations
fig_bpctr_stat = figure('position',[180,100,880,850]);

% ==================================================
% Left clicks, bpctr_opt = max
subplot(321)
if plot_type==0
    plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_left],'o-');
    hold on
    plot([data.bpctr_az.mean(freq_idx_wanted).max_left],'o-');
elseif plot_type==1
    errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_left],...
             [model.bpctr_az.std(model_freqI(freq_idx_wanted)).max_left],'o-');
    hold on
    errorbar([data.bpctr_az.mean(freq_idx_wanted).max_left],...
             [data.bpctr_az.std(freq_idx_wanted).max_left],'o-');
end
grid
legend('Model','Data','location','best');
title('Left, beam center = max')
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% Right clicks, bpctr_opt = max
subplot(322)
if plot_type==0
    plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_right],'o-');
    hold on
    plot([data.bpctr_az.mean(freq_idx_wanted).max_right],'o-');
elseif plot_type==1
    errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_right],...
             [model.bpctr_az.std(model_freqI(freq_idx_wanted)).max_right],'o-');
    hold on
    errorbar([data.bpctr_az.mean(freq_idx_wanted).max_right],...
             [data.bpctr_az.std(freq_idx_wanted).max_right],'o-');
end
grid
legend('Model','Data','location','best')
title('Right, beam center = max')
grid
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% ==================================================
% Left clicks, bpctr_opt = top
subplot(323)
if plot_type==0
    plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_left],'o-');
    hold on
    plot([data.bpctr_az.mean(freq_idx_wanted).top_left],'o-');
elseif plot_type==1
    errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_left],...
             [model.bpctr_az.std(model_freqI(freq_idx_wanted)).top_left],'o-');
    hold on
    errorbar([data.bpctr_az.mean(freq_idx_wanted).top_left],...
             [data.bpctr_az.std(freq_idx_wanted).top_left],'o-');
end
legend('Model','Data','location','best')
title('Left, beam center = top')
grid
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% Right clicks, bpctr_opt = top
subplot(324)
if plot_type==0
    plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_right],'o-');
    hold on
    plot([data.bpctr_az.mean(freq_idx_wanted).top_right],'o-');
elseif plot_type==1
    errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_right],...
             [model.bpctr_az.std(model_freqI(freq_idx_wanted)).top_right],'o-');
    hold on
    errorbar([data.bpctr_az.mean(freq_idx_wanted).top_right],...
             [data.bpctr_az.std(freq_idx_wanted).top_right],'o-');
end
legend('Model','Data','location','best')
title('Right, beam center = top')
grid
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% ==================================================
% Left clicks, bpctr_opt = ectr
subplot(325)
if plot_type==0
    plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_left],'o-');
    hold on
    plot([data.bpctr_az.mean(freq_idx_wanted).ectr_left],'o-');
elseif plot_type==1
    errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_left],...
             [model.bpctr_az.std(model_freqI(freq_idx_wanted)).ectr_left],'o-');
    hold on
    errorbar([data.bpctr_az.mean(freq_idx_wanted).ectr_left],...
             [data.bpctr_az.std(freq_idx_wanted).ectr_left],'o-');
end
legend('Model','Data','location','best')
title('Left, beam center = ectr')
grid
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% Right clicks, bpctr_opt = ectr
subplot(326)
if plot_type==0
    plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_right],'o-');
    hold on
    plot([data.bpctr_az.mean(freq_idx_wanted).ectr_right],'o-');
elseif plot_type==1
    errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_right],...
             [model.bpctr_az.std(model_freqI(freq_idx_wanted)).ectr_right],'o-');
    hold on
    errorbar([data.bpctr_az.mean(freq_idx_wanted).ectr_right],...
             [data.bpctr_az.std(freq_idx_wanted).ectr_right],'o-');
end
legend('Model','Data','location','best')
title('Right, beam center = ectr')
grid
set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
xlabel('Frequency (kHz)')
ylabel('Azimuth angle (deg)')

% Select ylim and set plot_type_str
for ii=1:6
    subplot(3,2,ii)
    if plot_type==0
        ylim([-15,15]);
        plot_type_str = 'meanonly';
    elseif plot_type==1
        ylim([-25,25]);
        plot_type_str = 'withstr';
    end
end


saveas(fig_bpctr_stat,fullfile(save_path,sprintf('%s_%s.fig',save_fname,plot_type_str)),'fig');
saveSameSize_res(fig_bpctr_stat,150,'file',...
                 fullfile(save_path,sprintf('%s_%s.png',save_fname,plot_type_str)),...
                 'format','png','renderer','painters');
epswrite(fullfile(save_path,sprintf('%s_%s.eps',save_fname,plot_type_str)));




%% MANOVA between data and model bpctr
% Check data for all [] --> make NaN
for iS=1:size(data.bpctr,1)
    for iF=1:size(data.bpctr,2)
        if isempty(data.bpctr(iS,iF).ectr_az)
            data.bpctr(iS,iF).ectr_az = NaN;
            data.bpctr(iS,iF).ectr_el = NaN;
        end
    end
end

% Check model for all [] --> make NaN
for iS=1:size(data.bpctr,1)
    for iF=1:size(data.bpctr,2)
        if isempty(model.bpctr(iS,iF).ectr_az)
            model.bpctr(iS,iF).ectr_az = NaN;
            model.bpctr(iS,iF).ectr_el = NaN;
        end
    end
end

% Data
for iF=1:length(data.param.freq_wanted)
    % mean of all beam center
    data_max(:,iF) = [data.bpctr(:,iF).max_az];
    data_top(:,iF) = [data.bpctr(:,iF).top_az];
    data_ectr(:,iF) = [data.bpctr(:,iF).ectr_az];
end

% Model
cnt = 0;
for iF=model_freqI
    cnt = cnt+1;
    model_max(:,cnt) = [model.bpctr(:,iF).max_az];
    model_top(:,cnt) = [model.bpctr(:,iF).top_az];
    model_ectr(:,cnt) = [model.bpctr(:,iF).ectr_az];
end

% MANOVA test
group = [zeros(300,1);ones(300,1)];
[d_max,p_max] = manova1([data_max;model_max],group);
[d_top,p_top] = manova1([data_top;model_top],group);
[d_ectr,p_ectr] = manova1([data_ectr;model_ectr],group);

[d_max,p_max] = manova1([data_max(:,3:4);model_max(:,3:4)],group);
[d_top,p_top] = manova1([data_top(:,3:4);model_top(:,3:4)],group);
[d_ectr,p_ectr] = manova1([data_ectr(:,3:4);model_ectr(:,3:4)],group);

