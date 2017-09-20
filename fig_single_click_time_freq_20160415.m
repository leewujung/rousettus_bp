% 2016 04 15  Plot time series and click spectrum of single channel

clear

ch_num = 33;  % mic #33 is GRAS 1/8" mic
freq_wanted = 35e3;
data_path = 'E:\fromDropbox\bp_processing\proc_output_rousettus_checked';
save_path = 'E:\fromDropbox\bp_processing\20160415_single_click_time_freq';

%% Select which click to use among all files and clicks
files = dir(fullfile(data_path,'*.mat'));
for iF=1:length(files)
    data_file = files(iF).name;
    ss = strsplit(data_file,'_');
    ss_p = strjoin(ss(3:4),'_');
    data = load(fullfile(data_path,data_file));
    call_num_total = size(data.proc.call_psd_dB_comp_re20uPa_withbp,1);
    for iC=1:call_num_total
        [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted,iC);
        [mm,mmidx] = max(call_dB);
        %         if call_dB(ch_num)>mm-1 && data.proc.chk_good_call(iC)==1
        if mmidx == ch_num && data.proc.chk_good_call(iC)==1
            fprintf('Click #%02d of file %s\n',iC,ss_p);
            data = load(fullfile(data_path,data_file));
            click = data.proc.call_align_short{iC,ch_num};
            t_click = (0:length(click)-1)/data.mic_data.fs;
            click_fft = fft(click);
            click_fft_dB = 20*log10(abs(click_fft));
            click_fft_dB_norm = click_fft_dB-max(click_fft_dB);
            freq_click = linspace(0,data.mic_data.fs,length(click_fft));
            
            ss_title = sprintf('bat %s, trial %s, click #%02d',ss{3},ss{4},iC);
            ss_file = sprintf('bat%s_trial%s_click%02d',ss{3},ss{4},iC);
            
            figure
            subplot(211)
            plot(t_click*1e6,click);
            xlabel('Time (us)');
            grid on
            title(ss_title);
            subplot(212)
            plot(freq_click/1e3,click_fft_dB_norm);
            grid on
            xlim([0 125]);
            ylim([-40 5])
            xlabel('Frequency (kHz)');
            saveas(gcf,fullfile(save_path,[ss_file,'.png']),'png');
        end
        
    end
end

%% Plot a particlar click: time series and spectrum
data_file = 'rousettus_20150825_34271_04_mic_data_bp_proc.mat';
data = load(fullfile(data_path,data_file));
call_num = 20;
click = data.proc.call_align_short{call_num,ch_num};
t_click = (0:length(click)-1)/data.mic_data.fs;
click_fft = fft(click);
click_fft_dB = 20*log10(abs(click_fft));
click_fft_dB_norm = click_fft_dB-max(click_fft_dB);
freq_click = linspace(0,data.mic_data.fs,length(click_fft));

psd_cmp = data.proc.call_psd_dB_comp_re20uPa_withbp{call_num,ch_num};
psd_cmp_norm = psd_cmp-max(psd_cmp);
freq_cmp = data.proc.call_freq_vec{call_num,ch_num};

ss = strsplit(data_file,'_');
ss_title = sprintf('bat %s, trial %s, click #%02d',ss{3},ss{4},call_num);

%% Time series and compensated spectrum (2 panel)
figure
subplot(211)
plot(t_click*1e6,click,'linewidth',2);
xlabel('Time (us)');
grid on
title(ss_title);
subplot(212)
plot(freq_cmp/1e3,psd_cmp_norm,'linewidth',2);
grid on
xlim([10 100]);
ylim([-30 5])
xlabel('Frequency (kHz)');
ylabel('Spectrum (dB)');

%% Time series, raw and compensated spectrum (3 panel)
figure
subplot(311)
plot(t_click*1e6,click,'linewidth',2);
xlabel('Time (us)');
grid on
title(ss_title);
subplot(312)
plot(freq_click/1e3,click_fft_dB_norm,'linewidth',2);
grid on
xlim([10 100]);
ylim([-30 5])
xlabel('Frequency (kHz)');
ylabel('Raw (dB)');
subplot(313)
plot(freq_cmp/1e3,psd_cmp_norm,'linewidth',2);
grid on
xlim([10 100]);
ylim([-30 5])
xlabel('Frequency (kHz)');
ylabel('Compensated (dB)');


%% Compare spectrum at increasingly off-axis angle from different channels
data_file = 'rousettus_20150825_34271_04_mic_data_bp_proc.mat';
data = load(fullfile(data_path,data_file));
call_num = 20;
ch_num = [33,30,28,23,19];

ss = strsplit(data_file,'_');
ss_title = sprintf('bat %s, trial %s, click #%02d',ss{3},ss{4},call_num);

fig_time = figure;
fig_freq = figure;
fig_freq_raw = figure;
for iCH=1:length(ch_num)
    click = data.proc.call_align_short{call_num,ch_num(iCH)};
    t_click = (0:length(click)-1)/data.mic_data.fs;

    psd_cmp = data.proc.call_psd_dB_comp_re20uPa_withbp{call_num,ch_num(iCH)};
    psd_cmp_norm = psd_cmp-max(psd_cmp);
    freq_cmp = data.proc.call_freq_vec{call_num,ch_num(iCH)};
    
    figure(fig_time)
    plot(t_click*1e6,click-0.1*(iCH-1),'linewidth',2);  % shift each channel by -0.1
    hold on
    
    figure(fig_freq_raw)
    plot(freq_cmp/1e3,psd_cmp,'linewidth',2);  % NO shifting across channels
    hold on

    figure(fig_freq)
    plot(freq_cmp/1e3,psd_cmp_norm-5*(iCH-1),'linewidth',2);  % shift each channel by -5
    hold on
end
figure(fig_time)
xlabel('Time (us)');
legend(num2str(ch_num'))
title(ss_title);

figure(fig_freq)
grid on
xlim([10 100]);
ylim([-50 5])
xlabel('Frequency (kHz)');
ylabel('Compensated & shifted (dB)');
legend(num2str(ch_num'))
title(ss_title);

figure(fig_freq_raw)
grid on
xlim([10 100]);
ylim([90 140])
xlabel('Frequency (kHz)');
ylabel('PSD (dB)');
legend(num2str(ch_num'))
title(ss_title);
