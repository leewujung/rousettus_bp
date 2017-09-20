% 2015 11 20  Effects of lowpass filtering clicks

clear
usrn = getenv('username');
% click#6 is better for channel 33 (on-axis)
ch34 = load(['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\sig_echo_test\rousettus_20150825_34271_02_click4_ch34.mat']);
ch26 = load(['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\sig_echo_test\rousettus_20150825_34271_02_click4_ch26.mat']);
load(['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\sig_echo_test\lpf_coef.mat']);

short_idx = 500:900;
ch34.call_short = ch34.call_align(short_idx);
w = tukeywin(length(ch34.call_short),0.25);

lpf_names = whos('lpf_*');

% sig_filt = nan(length(call_short),length(lpf_names));
% sm_env = sig_filt;
for iF=1:length(lpf_names)
    numfilt = eval(lpf_names(iF).name);
    sig_filt(:,iF) = filter(numfilt,1,ch34.call_align);
%     sig_filt(:,iF) = filtfilt(numfilt,1,ch34.call_align);
    sm_env(:,iF) = smooth(abs(hilbert(sig_filt(:,iF))),20);
end

sig_t = (0:size(sig_filt,1)-1)/ch34.fs;

figure
subplot(121)
plot(sig_t,bsxfun(@rdivide,sig_filt,max(sig_filt,[],1))+repmat(1:length(lpf_names),size(sig_filt,1),1));
subplot(122)
plot(sig_t,bsxfun(@rdivide,sm_env,max(sm_env,[],1))+repmat(1:length(lpf_names),size(sig_filt,1),1));

% Coherent summation enegry
ch34.call_fft = fft(ch34.call_short);
ch26.call_fft = fft(ch26.call_short);


% Percentage of energy included in short section



