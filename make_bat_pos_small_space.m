% 2015 11 12  Separate bat track if there is a gap in the middle

% Path/param
bat_pos_path = 'C:\Users\wlee76\Dropbox\0_ANALYSIS\bp_processing\bat_pos';
bat_pos_file = 'rousettus_20150825_34271_02_bat_pos.mat';

% Process track
bat = load(fullfile(bat_pos_path,bat_pos_file));
pos = cell2mat(bat.bat_pos);
pos = reshape(pos,length(pos),3,[]);
pos = nanmean(pos,3);  % raw bat track: mean of three points

notnanidx = find(~isnan(pos(:,1)));

% Find index for each track segment
jump_idx = find(diff(notnanidx)>10);
jump_idx = jump_idx(:)';
jump_idx = [0,jump_idx,length(notnanidx)];

if ~(jump_idx(1)==0 && jump_idx(2)==0)
    seg_idx = nan(length(jump_idx)-1,2);
    for iJ=1:length(jump_idx)-1
        seg_idx(iJ,:) = notnanidx([jump_idx(iJ)+1 jump_idx(iJ+1)]);
    end
    seg_idx(diff(seg_idx,1,2)<50,:) = [];
end

for iS=1:size(seg_idx,1)
    idx = ~ismember((1:length(pos)),seg_idx(iS,1):seg_idx(iS,2));
    bat_new = bat;
    for iM=1:3  % each marker
        bat_new.bat_pos{iM}(idx,:) = NaN;
    end
    save(fullfile(bat_pos_path,sprintf('%s_seg%02d.mat',strtok(bat_pos_file,'.'),iS)),'-struct','bat_new');
end
            

