% 2015 11 10  Make bat_pos file in the bp_proc format

% Path/param
usrn = getenv('username');
bat_pos_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing_large_room\bat_pos'];
mic_loc_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing_large_room\misc\mic_loc'];
pos_file = dir(fullfile(bat_pos_path,'*.csv'));
mic_wanted = [17,18,19,20,21,22,23];
sm_len = 25;
for iF=1:length(pos_file)
    pos = csvread(fullfile(bat_pos_path,pos_file(iF).name),1,0);
    pos = pos(:,1:3);  % only take first point
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
        seg_idx(diff(seg_idx,1,2)<100,:) = [];
        
        % Smooth track
        pos_sm(:,1) = smooth(pos(:,1),sm_len);
        pos_sm(:,2) = smooth(pos(:,2),sm_len);
        pos_sm(:,3) = smooth(pos(:,3),sm_len);
        
        % Estimate normal vector
        C = strsplit(pos_file(iF).name,'_');
        ddate = C{2};
        load(fullfile(mic_loc_path,[ddate,'_mic_loc.mat']));
        [U,S,V] = svd(mic_loc(mic_wanted,:));
        norm_vec = V(:,3);
        
        % Prepare to save
        ss = strsplit(pos_file(iF).name,'.csv');
        ss = ss{1};
        for iS=1:size(seg_idx,1)
            idx = seg_idx(iS,1):seg_idx(iS,2);
            A.bat_pos{1} = nan(size(pos));  % bat position
            A.bat_pos{1}(idx,:) = pos_sm(idx,:);
            A.bat_pos{2} = A.bat_pos{1};  % dummy columns
            A.bat_pos{3} = A.bat_pos{1};
            
            A.nn = nan(size(pos));  % bat head normal
            A.nn(idx,:) = repmat(norm_vec',length(idx),1);
            
            hv = diff(A.bat_pos{1});  % estimated head aim vec
            hv(idx(end),:) = hv(idx(end-1),:);  % use the last one to substitute
            hv = hv./repmat(sqrt(diag(hv*hv')),1,3);  % normalize estimated head aim vec
            
            A.head_vec = nan(size(pos));  % bat head aim
            A.head_vec(idx,:) = hv(idx,:);
            save(fullfile(bat_pos_path,[ss,sprintf('_seg%02d.mat',iS)]),'-struct','A');
        end
        
    end
    
    clearvars -except bat_pos_path mic_loc_path pos_file mic_wanted sm_len
    
end

