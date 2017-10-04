% 2017 09 25  Compare the beam center location for data and model, use all
%             three possible defitions of beam center
% 2017 09 26  MANOVA to see if model and data bpctr can be consider "the
%             same" as a whole across multiple frequencies
% 2017 10 04  -- Update to include piston model results
%             -- Do stat test on bpctr location at each frequency individually

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
piston_path = 'model_piston_multifreq_indiv_bpctr_cntr_20171003_std1.0';
piston_file = 'model_piston_multifreq_indiv_bpctr_cntr_20171003_std1.0_results.mat';

model = load(fullfile(data_base_path,results_path,model_path,model_file));
data = load(fullfile(data_base_path,results_path,data_path,data_file));
piston = load(fullfile(data_base_path,results_path,piston_path,piston_file));

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
plot_type = 1; % 0-'meanonly', 1-'withstd'



%% Plot mean beam center az locations
fig_bpctr_stat = figure('position',[180,100,880,850]);

for ii=1:3
    % Left clicks, bpctr_opt = max
    subplot(3,2,(ii-1)*2+1)
    if plot_type==0
        if ii==1  % bpctr_opt = max
            hd  = plot([data.bpctr_az.mean(freq_idx_wanted).max_left],'o-');
            hold on
            hph = plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_left],'o-');
            hpi = plot([piston.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_left],'o-');
        elseif ii==2  % bpctr_opt = top
            hd  = plot([data.bpctr_az.mean(freq_idx_wanted).top_left],'o-');
            hold on
            hph = plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_left],'o-');
            hpi = plot([piston.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_left],'o-');
        elseif ii==3  % bpctr_opt = ectr
            hd  = plot([data.bpctr_az.mean(freq_idx_wanted).ectr_left],'o-');
            hold on
            hph = plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_left],'o-');
            hpi = plot([piston.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_left],'o-');
        end
    elseif plot_type==1
        if ii==1  % bpctr_opt = max
            hd  = errorbar([data.bpctr_az.mean(freq_idx_wanted).max_left],...
                           [data.bpctr_az.std(freq_idx_wanted).max_left],'o-');
            hold on 
            hph = errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_left],...
                           [model.bpctr_az.std(model_freqI(freq_idx_wanted)).max_left],'o-');
            hpi = errorbar([piston.bpctr_az.mean(freq_idx_wanted).max_left],...
                           [piston.bpctr_az.std(freq_idx_wanted).max_left],'o-');
        elseif ii==2  % bpctr_opt = top
            hd  = errorbar([data.bpctr_az.mean(freq_idx_wanted).top_left],...
                           [data.bpctr_az.std(freq_idx_wanted).top_left],'o-');
            hold on 
            hph = errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_left],...
                           [model.bpctr_az.std(model_freqI(freq_idx_wanted)).top_left],'o-');
            hpi = errorbar([piston.bpctr_az.mean(freq_idx_wanted).top_left],...
                           [piston.bpctr_az.std(freq_idx_wanted).top_left],'o-');
        elseif ii==3  % bpctr_opt = ectr
            hd  = errorbar([data.bpctr_az.mean(freq_idx_wanted).ectr_left],...
                           [data.bpctr_az.std(freq_idx_wanted).ectr_left],'o-');
            hold on 
            hph = errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_left],...
                           [model.bpctr_az.std(model_freqI(freq_idx_wanted)).ectr_left],'o-');
            hpi = errorbar([piston.bpctr_az.mean(freq_idx_wanted).ectr_left],...
                           [piston.bpctr_az.std(freq_idx_wanted).ectr_left],'o-');
        end
    end
    grid
    legend([hd,hpi,hph],{'Data','Piston','Phased',''},'location','best');
    title('Left, beam center = max')
    set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
    xlabel('Frequency (kHz)')
    ylabel('Azimuth angle (deg)')

    % Right clicks, bpctr_opt = max
    subplot(3,2,(ii-1)*2+2)
    if plot_type==0
        if ii==1  % bpctr_opt = max
            hd  = plot([data.bpctr_az.mean(freq_idx_wanted).max_right],'o-');
            hold on
            hph = plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_right],'o-');
            hpi = plot([piston.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_right],'o-');
        elseif ii==2  % bpctr_opt = top
            hd  = plot([data.bpctr_az.mean(freq_idx_wanted).top_right],'o-');
            hold on
            hph = plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_right],'o-');
            hpi = plot([piston.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_right],'o-');
        elseif ii==3  % bpctr_opt = ectr
            hd  = plot([data.bpctr_az.mean(freq_idx_wanted).ectr_right],'o-');
            hold on
            hph = plot([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_right],'o-');
            hpi = plot([piston.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_right],'o-');
        end
    elseif plot_type==1
        if ii==1  % bpctr_opt = max
            hd  = errorbar([data.bpctr_az.mean(freq_idx_wanted).max_right],...
                           [data.bpctr_az.std(freq_idx_wanted).max_right],'o-');
            hold on 
            hph = errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).max_right],...
                           [model.bpctr_az.std(model_freqI(freq_idx_wanted)).max_right],'o-');
            hpi = errorbar([piston.bpctr_az.mean(freq_idx_wanted).max_right],...
                           [piston.bpctr_az.std(freq_idx_wanted).max_right],'o-');
        elseif ii==2  % bpctr_opt = top
            hd  = errorbar([data.bpctr_az.mean(freq_idx_wanted).top_right],...
                           [data.bpctr_az.std(freq_idx_wanted).top_right],'o-');
            hold on 
            hph = errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).top_right],...
                           [model.bpctr_az.std(model_freqI(freq_idx_wanted)).top_right],'o-');
            hpi = errorbar([piston.bpctr_az.mean(freq_idx_wanted).top_right],...
                           [piston.bpctr_az.std(freq_idx_wanted).top_right],'o-');
        elseif ii==3  % bpctr_opt = ectr
            hd  = errorbar([data.bpctr_az.mean(freq_idx_wanted).ectr_right],...
                           [data.bpctr_az.std(freq_idx_wanted).ectr_right],'o-');
            hold on 
            hph = errorbar([model.bpctr_az.mean(model_freqI(freq_idx_wanted)).ectr_right],...
                           [model.bpctr_az.std(model_freqI(freq_idx_wanted)).ectr_right],'o-');
            hpi = errorbar([piston.bpctr_az.mean(freq_idx_wanted).ectr_right],...
                           [piston.bpctr_az.std(freq_idx_wanted).ectr_right],'o-');
        end
    end
    grid
    legend([hd,hpi,hph],{'Data','Piston','Phased',''},'location','best');
    title('Left, beam center = max')
    set(gca,'xtick',1:9,'xticklabel',[model.freq.all(model_freqI(freq_idx_wanted))'/1e3])
    xlabel('Frequency (kHz)')
    ylabel('Azimuth angle (deg)')

end


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




%% Paired KS 2-sample test for bpctr at each frequency
for iF = 1:length(freq_idx_wanted)
    % KS 2-sample test
    % bpctr_opt = max
    [ks.h_data_piston(iF,1),ks.p_data_piston(iF,1)] = ...
        kstest2([data.bpctr(:,freq_idx_wanted(iF)).max_az],...
                [piston.bpctr(:,model_freqI(freq_idx_wanted(iF))).max_az]);
    [ks.h_data_phase(iF,1),ks.p_data_phase(iF,1)] = ...
        kstest2([data.bpctr(:,freq_idx_wanted(iF)).max_az],...
                [model.bpctr(:,model_freqI(freq_idx_wanted(iF))).max_az]);
    % bpctr_opt = top
    [ks.h_data_piston(iF,2),ks.p_data_piston(iF,2)] = ...
        kstest2([data.bpctr(:,freq_idx_wanted(iF)).top_az],...
                [piston.bpctr(:,model_freqI(freq_idx_wanted(iF))).top_az]);
    [ks.h_data_phase(iF,2),ks.p_data_phase(iF,2)] = ...
        kstest2([data.bpctr(:,freq_idx_wanted(iF)).top_az],...
                [model.bpctr(:,model_freqI(freq_idx_wanted(iF))).top_az]);
    % bpctr_opt = ectr
    [ks.h_data_piston(iF,3),ks.p_data_piston(iF,3)] = ...
        kstest2([data.bpctr(:,freq_idx_wanted(iF)).ectr_az],...
                [piston.bpctr(:,model_freqI(freq_idx_wanted(iF))).ectr_az]);
    [ks.h_data_phase(iF,3),ks.p_data_phase(iF,3)] = ...
        kstest2([data.bpctr(:,freq_idx_wanted(iF)).ectr_az],...
                [model.bpctr(:,model_freqI(freq_idx_wanted(iF))).ectr_az]);

    % Mann-Whitney U-test
    % bpctr_opt = max
    [mw.p_data_piston(iF,1),mw.h_data_piston(iF,1),mw.stat_data_piston(iF,1)] = ...
        ranksum([data.bpctr(:,freq_idx_wanted(iF)).max_az],...
                [piston.bpctr(:,model_freqI(freq_idx_wanted(iF))).max_az]);
    [mw.p_data_phase(iF,1),mw.h_data_phase(iF,1),mw.stat_data_phase(iF,1)] = ...
        ranksum([data.bpctr(:,freq_idx_wanted(iF)).max_az],...
                [model.bpctr(:,model_freqI(freq_idx_wanted(iF))).max_az]);
    % bpctr_opt = top
    [mw.p_data_piston(iF,2),mw.h_data_piston(iF,2),mw.stat_data_piston(iF,2)] = ...
        ranksum([data.bpctr(:,freq_idx_wanted(iF)).top_az],...
                [piston.bpctr(:,model_freqI(freq_idx_wanted(iF))).top_az]);
    [mw.p_data_phase(iF,2),mw.h_data_phase(iF,2),mw.stat_data_phase(iF,2)] = ...
        ranksum([data.bpctr(:,freq_idx_wanted(iF)).top_az],...
                [model.bpctr(:,model_freqI(freq_idx_wanted(iF))).top_az]);
    % bpctr_opt = ectr
    [mw.p_data_piston(iF,3),mw.h_data_piston(iF,3),mw.stat_data_piston(iF,3)] = ...
        ranksum([data.bpctr(:,freq_idx_wanted(iF)).ectr_az],...
                [piston.bpctr(:,model_freqI(freq_idx_wanted(iF))).ectr_az]);
    [mw.p_data_phase(iF,3),mw.h_data_phase(iF,3),mw.stat_data_phase(iF,3)] = ...
        ranksum([data.bpctr(:,freq_idx_wanted(iF)).ectr_az],...
                [model.bpctr(:,model_freqI(freq_idx_wanted(iF))).ectr_az]);
end




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

