% 2015 11 10  Extract channels from large room mic_data

% Param
seq_path = 'C:\Users\wlee76\Dropbox\0_ANALYSIS\bp_processing_large_room\misc\mic_seq';
mic_data_path = 'C:\Users\wlee76\Dropbox\0_ANALYSIS\bp_processing_large_room\mic_data';
mic_detect_path = 'C:\Users\wlee76\Dropbox\0_ANALYSIS\bp_processing_large_room\mic_detect';
ddate = {'20141204','20141205','20141219'};

for iD = 1:length(ddate)
    match = csvread(fullfile(seq_path,[ddate{iD},'_mic_seq.csv']),1,0);
    match = sort(match(:,2));
    % mic_data files
    ff_data = dir(fullfile(mic_data_path,['rousettus_',ddate{iD},'*.mat']));
    for iF = 1:length(ff_data)
        A = load(fullfile(mic_data_path,ff_data(iF).name));
        A.sig = A.sig(:,match);
        save(fullfile(mic_data_path,ff_data(iF).name),'-struct','A');
    end
    % _detect files
    ff_detect = dir(fullfile(mic_detect_path,['rousettus_',ddate{iD},'*.mat']));
    for iF = 1:length(ff_detect)
        B = load(fullfile(mic_detect_path,ff_detect(iF).name));
        B.sig_rough = B.sig_rough(:,match);
        B.num_ch_in_file = length(match);
        save(fullfile(mic_detect_path,ff_detect(iF).name),'-struct','B');
    end
end