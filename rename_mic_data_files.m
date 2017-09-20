% 2015 11 08  Rename mic_data, mic_data_detect files

% % To add date to the filename
% base_path = 'E:\0_WJLEE_DATA\Beampattern Rousettus';
% ddate = {'20141204','20141205','20141219'};
% 
% for iD=1:length(ddate)
%     rel_path = sprintf('%s_rousettus\\%s_rousettus_mic_matfile',ddate{iD},ddate{iD});
%     ff = dir(fullfile(base_path,rel_path,'*.mat'));
%     for iF=1:length(ff)
%         ss = ff(iF).name;
%         C = strsplit(ss,'_');
%         C{1} = sprintf('%s_%s','rousettus',ddate{iD});
%         C{2} = sprintf('%02d',str2double(C{2}));
%         ss_new = strjoin(C,'_');
%         movefile(fullfile(base_path,rel_path,ss),...
%                  fullfile(base_path,rel_path,ss_new));
%     end
% end


% Add bat name into filename
base_path = 'C:\Users\Wu-Jung Lee\Dropbox\0_ANALYSIS\bp_processing_large_room\mic_data';
ddate = {'20141204','20141205','20141219'};
indiv{1} = {(1:7),(8:16),(17:19)};
bat{1} = {'BK','GR','LU'};
indiv{2} = {(1:7),(8:10)};
bat{2} = {'BK','GR'};
indiv{3} = {(2:5),(6:9),(10:13),(14:19)};
bat{3} = {'RD','BK','GR','BL'};

for iD=1:length(ddate)
    ff = dir(fullfile(base_path,['rousettus_',ddate{iD},'*.mat']));
    for iF=1:length(ff)
        ss = ff(iF).name;
        C = strsplit(ss,'_');
        CS = C;
        C{1} = CS{2};
        C{2} = CS{1};
        num = str2double(C{3});
        C(4+1:length(C)+1) = deal(C(4:end));
        for iB=1:length(indiv{iD})
            if ismember(num,indiv{iD}{iB})
                C{4} = bat{iD}{iB};
            end
        end
        ss_new = strjoin(C,'_');
        movefile(fullfile(base_path,ss),fullfile(base_path,ss_new));
    end
end




