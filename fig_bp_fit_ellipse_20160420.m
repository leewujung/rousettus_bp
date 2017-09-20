% 2016 04 20  Plot single click 35 kHz beam pattern before and after
%             shift/tilt rotation and also with the best-fitting ellipse


clear
usrn = getenv('username');
if strcmp(usrn,'Wu-Jung')
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\export_fig-master');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\export_fig-master']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

% Set up various paths
[~,script_name,~] = fileparts(mfilename('fullpath'));
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_processing';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing'];
end

save_root = 'analysis_results_figs';
save_path = fullfile(base_path,save_root,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

data_path = 'rotate_all_click_20160419';
data_file = 'rotate_all_click_20160419_all_click_rotated_36134.mat';
D = load(fullfile(base_path,save_root,data_path,data_file));

% Select file & click
bat = '36134';
trial = 2;
click = 20;
fname = sprintf('rousettus_20150825_%s_%02d_mic_data_bp_proc.mat',bat,trial);
cvec = 0:-3:-30;
az_plot_range = 120;

% Find corresponding indices
for iB=1:length(D.bp_processed_file_all)
    if strcmp(D.bp_processed_file_all(iB).name,fname);
        break
    end
end
[~,iC] = ismember(click,D.raw_meas.good_call_idx{iB});


fig = figure;
set(fig,'position',[140 100 1100 480])

V = D.rotate_data{iB}(iC).raw;
V.vq_norm(abs(V.azq)>=az_plot_range) = NaN;

subplot(121)  % raw
axesm('eckert4');
axis off
% contourf(V.xq,V.yq,V.vq_norm,cvec(2:end),'w','linewidth',1);  % this is equivalent to using contourfm & azq/elq
contour(V.xq,V.yq,V.vq_norm,cvec(2:end),'fill','on');  % this is equivalent to using contourfm & azq/elq
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',az_plot_range*[-1 1]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
colorbar('southoutside','ticks',fliplr(cvec));
colormap(parula(length(cvec)-1));
caxis([cvec(end), 0])
tightmap
title('Raw')

V = D.rotate_data{iB}(iC).rot_elpctr_tilt;
V.vq_norm(V.outofbnd_azel_idx_q) = NaN;
V.vq_norm(abs(V.azq)>=az_plot_range) = NaN;
E = V.E;
c3db_xy = E.c3db_xy;
xmin = min(c3db_xy(:,1))-range(c3db_xy(:,1)*0.5);
xmax = max(c3db_xy(:,1))+range(c3db_xy(:,1)*0.5);
ymin = min(c3db_xy(:,2))-range(c3db_xy(:,1)*0.5);
ymax = max(c3db_xy(:,2))+range(c3db_xy(:,1)*0.5);

subplot(122)  % shift/tilt final
axesm('eckert4');
axis off
% contourf(V.xq,V.yq,V.vq_norm,cvec(2:end),'w','linewidth',1);  % this is equivalent to using contourfm & azq/elq
contour(V.xq,V.yq,V.vq_norm,cvec(2:end),'fill','on');  % this is equivalent to using contourfm & azq/elq
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',az_plot_range*[-1 1]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-','glinewidth',1);
hold on
fit_df = ezplot(E.eqt,[xmin,xmax,ymin,ymax]);  % plot ellipse
set(fit_df,'linecolor','b','linewidth',2,'linestyle','-');
plot(E.x0,E.y0,'r*','linewidth',1.5,'markersize',10);
colorbar('southoutside','ticks',fliplr(cvec));
colormap(parula(length(cvec)-1));
caxis([cvec(end), 0])
tightmap
title('Shift and tilted')

set(fig,'color','w');

bat = '36134';
trial = 2;
click = 20;

ss_file = sprintf('%s_%s_%02d_c%02d_az%03ddeg',...
    script_name,bat,trial,click,az_plot_range);
saveas(fig,fullfile(save_path,[ss_file,'.fig']),'fig');
% export_fig(fullfile(save_path,ss_file),'-q101','-eps');  % high quality EPS
saveSameSize(fig,'file',fullfile(save_path,ss_file),...
    'format','png','renderer','painters');




