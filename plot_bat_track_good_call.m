% 2015 12 04  Plot all bat tracks and good call location

clear
usrn = getenv('username');

% Bat data path
base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
bat_proc_dir = './proc_output';
bat_proc_file = dir(fullfile(base_path,bat_proc_dir,'rousettus_20150825_*.mat'));
freq = 35e3;

% Load bat track and call location
A = load(fullfile(base_path,bat_proc_dir,bat_proc_file(2).name));
track = nan(length(bat_proc_file),size(A.track.track_interp,1),size(A.track.track_interp,2));
call_loc = cell(length(bat_proc_file),1);
call_loc_idx_on_track = cell(length(bat_proc_file),1);
call_max_ch_idx = cell(length(bat_proc_file),1);
call_max_azel = cell(length(bat_proc_file),1);
call_all_azel = cell(length(bat_proc_file),1);
for iB = 1:length(bat_proc_file)
    A = load(fullfile(base_path,bat_proc_dir,bat_proc_file(iB).name));
    track(iB,:,:) = A.track.track_interp;
    good_call_idx = find(A.proc.chk_good_call==1);
    call_loc{iB} = A.proc.bat_loc_at_call(good_call_idx,:);
    call_loc_idx_on_track{iB} = A.proc.call_loc_on_track_interp(good_call_idx);
end
mic_loc = A.mic_loc;
call_loc = cell2mat(call_loc);  % collapse all call_loc together

% Mic numbers for each sides
left_idx = [1,16,6,3,2,8];
top_idx = [13,15,9,11,14,12,10];
bottom_idx = [27,4,26,5,17,24,7];
right_idx = [18,21,23,22,19,20];
front_idx = [29,30,28,34,33,25,32,31];

% Plot good call location
fig_track = figure;
corder = get(gca,'colororder');
hmic_top = plot3(mic_loc(top_idx,1),mic_loc(top_idx,2),mic_loc(top_idx,3),'ko','markerfacecolor','k');
hold on
hmic_bottom = plot3(mic_loc(bottom_idx,1),mic_loc(bottom_idx,2),mic_loc(bottom_idx,3),'ko','markerfacecolor','k');
hmic_left = plot3(mic_loc(left_idx,1),mic_loc(left_idx,2),mic_loc(left_idx,3),'ko');
hmic_right = plot3(mic_loc(right_idx,1),mic_loc(right_idx,2),mic_loc(right_idx,3),'ko');
hmic_front = plot3(mic_loc(front_idx,1),mic_loc(front_idx,2),mic_loc(front_idx,3),'ko');
grid on
axis equal
hold on
for iB=1:size(track,1)
    tr = squeeze(track(iB,:,:));
    htrack = plot3(tr(:,1),tr(:,2),tr(:,3),'color',corder(1,:));
end
% hcall = plot3(call_loc(:,1),call_loc(:,2),call_loc(:,3),'.','color',corder(2,:));
hcall = plot3(call_loc(:,1),call_loc(:,2),call_loc(:,3),'.','color','r');
legend([hmic_top,hmic_left,htrack,hcall],{'Mic (ceiling/floor)','Mic (left/right/front)','Bat track','Good calls'},...
       'location','eastoutside');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');

% Plot good call location in X-Y view
fig_track = figure;
corder = get(gca,'colororder');
hmic_left = plot3(mic_loc(left_idx,1),mic_loc(left_idx,2),mic_loc(left_idx,3),'ko');
hold on
hmic_right = plot3(mic_loc(right_idx,1),mic_loc(right_idx,2),mic_loc(right_idx,3),'ko');
hmic_front = plot3(mic_loc(front_idx,1),mic_loc(front_idx,2),mic_loc(front_idx,3),'ko');
grid on
axis equal
hold on
for iB=1:size(track,1)
    tr = squeeze(track(iB,:,:));
    htrack = plot3(tr(:,1),tr(:,2),tr(:,3),'color',corder(1,:));
end
% hcall = plot3(call_loc(:,1),call_loc(:,2),call_loc(:,3),'.','color',corder(2,:));
hcall = plot3(call_loc(:,1),call_loc(:,2),call_loc(:,3),'.','color','r');
legend([hmic_left,htrack,hcall],{'Mic (left/right/front)','Bat track','Good calls'},...
       'location','eastoutside');
view([0 0 1]);
xlabel('X (m)'); ylabel('Y (m)');

