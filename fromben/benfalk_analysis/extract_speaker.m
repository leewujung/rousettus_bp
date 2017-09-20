clear;

date_data='20150826';

trials=1:11;
speaker_pos=cell(length(trials),1);
for npos=1:length(trials)
  fn=['F:\small_space_beampattern\Test\' date_data...
    '\Session 1\Trial' num2str(trials(npos),'%2.2d') '.c3d'];
  [point_array, frame_rate, trig_rt, trig_sig, start_f, end_f] = lc3d( fn );
  point_names=cellfun(@(c) c.name,point_array,'uniformoutput',0);
  
  idx=cellfun('isempty',strfind(point_names,'*'));
  speaker_pos{npos}=cell2mat(...
    cellfun(@(c) c.traj(1,:)/1e3,point_array(idx),'uniformoutput',0)...
    );
end

figure(2); clf;
for k=1:length(trials)
  hold on;
  scatter3(speaker_pos{k}(:,1),speaker_pos{k}(:,2),speaker_pos{k}(:,3),'.')
end
axis equal; grid on;

save(['speaker_pos_' date_data '.mat'],'speaker_pos');