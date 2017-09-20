% 2015 11 11  Mkae mic_loc file for large room rousettus recording

clear

base_path = 'C:\Users\wlee76\Dropbox\0_ANALYSIS\bp_processing_large_room\misc';
mic_mark_path = 'mic_marking_stuff';
mic_loc_path = 'mic_loc';
mic_vec_path = 'mic_vec';
mic_seq_path = 'mic_seq';

% date: 20141031
dd = '20141031';
pos1 = csvread(fullfile(base_path,mic_mark_path,'20141031_mic_loc_set1_xyzpts.csv'),1,0);
pos1 = reshape(pos1(1,:),3,[])';
pos2 = csvread(fullfile(base_path,mic_mark_path,'20141031_mic_loc_set2_xyzpts.csv'),1,0);
pos2 = reshape(pos2(1,:),3,[])';
pos = [pos1;pos2];
vec1 = csvread(fullfile(base_path,mic_mark_path,'20141031_mic_vec_set1_xyzpts.csv'),1,0);
vec1 = reshape(vec1(1,:),3,[])';
vec2 = csvread(fullfile(base_path,mic_mark_path,'20141031_mic_vec_set2_xyzpts.csv'),1,0);
vec2 = reshape(vec2(1,:),3,[])';
vec = [vec1;vec2];
vec = vec(2:2:end,:)-vec(1:2:end,:);  % front-base for mic_vec

mm = csvread(fullfile(base_path,mic_seq_path,'20141031_mic_seq.csv'),1,0);
[mm,~] = sortrows(mm,2);  % sort according to channel sequence

mic_loc = nan(32,3);
mic_loc(mm(:,2),:) = pos(mm(:,1),:);  % fill in mic_loc at corresponding channel
mic_vec = nan(32,3);
mic_vec(mm(:,2),:) = vec(mm(:,1),:);  % fill in mic_loc at corresponding channel

save(fullfile(base_path,mic_loc_path,[dd,'_mic_loc.mat']),'mic_loc');
save(fullfile(base_path,mic_vec_path,[dd,'_mic_vec.mat']),'mic_vec');

figure
plot3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),'o');
hold on
text(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),num2str([1:32]'));
quiver3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),mic_vec(:,3),mic_vec(:,1),mic_vec(:,2));
title(num2str(dd));
grid on
axis equal
saveas(gcf,fullfile(base_path,mic_loc_path,[dd,'_mic_loc_vec.fig']));

clearvars -except base_path mic_*_path

% date: 20141204
dd = '20141204';
pos = csvread(fullfile(base_path,mic_mark_path,'20141204_mic_loc_xyzpts.csv'),1,0);
pos = reshape(pos(1,1:90),3,[])';
vec = csvread(fullfile(base_path,mic_mark_path,'20141204_mic_vec_xyzpts.csv'),1,0);
vec = reshape(vec(1,:),3,[])';
vec = vec(2:2:end,:)-vec(1:2:end,:);  % front-base for mic_vec

mm = csvread(fullfile(base_path,mic_seq_path,'20141204_mic_seq.csv'),1,0);
[mm,~] = sortrows(mm,2);  % sort according to channel sequence

mic_loc = nan(32,3);
mic_loc(mm(:,2),:) = pos(mm(:,1),:);  % fill in mic_loc at corresponding channel
mic_vec = nan(32,3);
mic_vec(mm(:,2),:) = vec(mm(:,1),:);  % fill in mic_loc at corresponding channel

save(fullfile(base_path,mic_loc_path,[dd,'_mic_loc.mat']),'mic_loc');
save(fullfile(base_path,mic_vec_path,[dd,'_mic_vec.mat']),'mic_vec');

figure
plot3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),'o');
hold on
text(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),num2str([1:32]'));
quiver3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),mic_vec(:,3),mic_vec(:,1),mic_vec(:,2));
title(num2str(dd));
grid on
axis equal
saveas(gcf,fullfile(base_path,mic_loc_path,[dd,'_mic_loc_vec.fig']));

clearvars -except base_path mic_*_path

% date: 20141205
dd = '20141205';
pos = csvread(fullfile(base_path,mic_mark_path,'20141205_mic_loc_xyzpts.csv'),1,0);
pos_new = pos(1,:);
pos_new(isnan(pos(1,:))) = pos(3,isnan(pos(1,:)));
pos = reshape(pos_new,3,[])';
vec = csvread(fullfile(base_path,mic_mark_path,'20141205_mic_vec_xyzpts.csv'),1,0);
vec_new = vec(1,:);
vec_new(isnan(vec(1,:))) = vec(3,~isnan(vec(3,:)));
vec = reshape(vec_new(1,:),3,[])';
vec = vec(2:2:end,:)-vec(1:2:end,:);  % front-base for mic_vec

mm_loc = csvread(fullfile(base_path,mic_seq_path,'20141205_mic_seq.csv'),1,0);
[mm_loc,~] = sortrows(mm_loc,2);  % sort according to channel sequence
mm_vec = csvread(fullfile(base_path,mic_seq_path,'20141205_mic_seq_vec.csv'),1,0);
[mm_vec,~] = sortrows(mm_vec,2);  % sort according to channel sequence

mic_loc = nan(32,3);
mic_loc(mm_loc(:,2),:) = pos(mm_loc(:,1),:);  % fill in mic_loc at corresponding channel
mic_vec = nan(32,3);
mic_vec(mm_vec(:,2),:) = vec(mm_vec(:,1),:);  % fill in mic_loc at corresponding channel

save(fullfile(base_path,mic_loc_path,[dd,'_mic_loc.mat']),'mic_loc');
save(fullfile(base_path,mic_vec_path,[dd,'_mic_vec.mat']),'mic_vec');

figure
plot3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),'o');
hold on
text(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),num2str([1:32]'));
quiver3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),mic_vec(:,3),mic_vec(:,1),mic_vec(:,2));
title(num2str(dd));
grid on
axis equal
saveas(gcf,fullfile(base_path,mic_loc_path,[dd,'_mic_loc_vec.fig']));

clearvars -except base_path mic_*_path

% date: 20141219
dd = '20141219';
pos = csvread(fullfile(base_path,mic_mark_path,'20141219_mic_loc_xyzpts.csv'),1,0);
pos = reshape(pos(1,:),3,[])';
vec = csvread(fullfile(base_path,mic_mark_path,'20141219_mic_vec_xyzpts.csv'),1,0);
vec = reshape(vec(1,:),3,[])';
vec = vec(2:2:end,:)-vec(1:2:end,:);  % front-base for mic_vec

mm_loc = csvread(fullfile(base_path,mic_seq_path,'20141219_mic_seq.csv'),1,0);
[mm_loc,~] = sortrows(mm_loc,2);  % sort according to channel sequence
mm_vec = csvread(fullfile(base_path,mic_seq_path,'20141219_mic_seq_vec.csv'),1,0);
[mm_vec,~] = sortrows(mm_vec,2);  % sort according to channel sequence

mic_loc = nan(32,3);
mic_loc(mm_loc(:,2),:) = pos(mm_loc(:,1),:);  % fill in mic_loc at corresponding channel
mic_vec = nan(32,3);
mic_vec(mm_vec(:,2),:) = vec(mm_vec(:,1),:);  % fill in mic_loc at corresponding channel

save(fullfile(base_path,mic_loc_path,[dd,'_mic_loc.mat']),'mic_loc');
save(fullfile(base_path,mic_vec_path,[dd,'_mic_vec.mat']),'mic_vec');

figure
plot3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),'o');
hold on
text(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),num2str([1:32]'));
quiver3(mic_loc(:,3),mic_loc(:,1),mic_loc(:,2),mic_vec(:,3),mic_vec(:,1),mic_vec(:,2));
title(num2str(dd));
grid on
axis equal
saveas(gcf,fullfile(base_path,mic_loc_path,[dd,'_mic_loc_vec.fig']));

clearvars -except base_path mic_*_path

