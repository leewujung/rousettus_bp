% extract_audio
% 
% 
% 
% bat 36134 
% 
% 2015.08.25
% 2 3 
% 
% 2015.09.10 
% 2 3
clear;

frm_rate=15;


file_dir='..\mic_recordings\20150825_calib_mic\rousettus_36134_matfile\';
% file_dir='..\mic_recordings\20150910_calib_mic\rousettus_36134_matfile\';

filename='rousettus_36134_2_mic_data.mat';
% filename='rousettus_36134_3_mic_data.mat';
% filename='rousettus_361342_mic_data.mat';
% filename='rousettus_361343_mic_data.mat';

load([file_dir filename]);

[~,ch]=max(max(sig));

[b,a]=butter(6,90e3/(fs/2),'low');

sig_filt=filtfilt(b,a,sig(:,ch));

audiowrite([file_dir filename(1:end-3) 'wav'],sig_filt,fs*frm_rate/200);