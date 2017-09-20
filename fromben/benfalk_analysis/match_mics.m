clear

fn='F:\small_space_beampattern\Test\20150819\Session 1\Trial02.c3d';

[point_array, frame_rate, trig_rt, trig_sig, start_f, end_f] = lc3d( fn );
all_points=cellfun(@(c) c.traj(1,:) , point_array,'uniformoutput',0);

pts=cell2mat(all_points)./1e3;


MP=load('mic_pos_20150820_S1_T1.mat');
early_MP=MP.mic_pos;


figure(1); clf; hold on;
scatter3(early_MP(:,1),early_MP(:,2),early_MP(:,3),'o')
grid on; axis equal;

scatter3(pts(:,1),pts(:,2),pts(:,3),'o')


for pp=1:length(pts)
  D=distance(early_MP,pts(pp,:));
  [MD(pp),ii(pp)]=min(D);
end

[mdsort,isort]=sort(MD);
iisort=ii(isort);

mic_pts=pts(iisort(1:3*32),:);


figure(2); clf; hold on;
scatter3(mic_pts(:,1),mic_pts(:,2),mic_pts(:,3),'o')
grid on; axis equal;

scatter3(early_MP(:,1),early_MP(:,2),early_MP(:,3),'o')

