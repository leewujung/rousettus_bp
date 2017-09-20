% 2015 11 08  Rename bat_pos files from large room recording


% % Add date into filename
% base_path = 'E:\0_WJLEE_DATA\Beampattern Rousettus';
% ddate = {'20141204','20141205','20141219'};
% 
% for iD=1:length(ddate)
%     rel_path = sprintf('%s_rousettus\\analyzed video files',ddate{iD});
%     ff = dir(fullfile(base_path,rel_path,'*xyzpts.csv'));
%     for iF=1:length(ff)
%         ss = ff(iF).name;
%         if iD==1
%             ss_new = sprintf('%s_rousettus_%02d_bat_pos.csv',ddate{iD},iF+1);
%         elseif iD==2
%             ss_new = sprintf('%s_rousettus_%02d_bat_pos.csv',ddate{iD},iF);
%         elseif iD==3
%             if iF<14
%                 ss_new = sprintf('%s_rousettus_%02d_bat_pos.csv',ddate{iD},iF+1);
%             else
%                 ss_new = sprintf('%s_rousettus_%02d_bat_pos.csv',ddate{iD},iF+2);
%             end
%         end
%         copyfile(fullfile(base_path,rel_path,ss),...
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