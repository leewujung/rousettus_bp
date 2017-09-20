plot_mic_nums=1;


mic_proc_dir='..\proc_output\';
fname='rousettus_20150825_34271_02_mic_data_bp_proc.mat';

bpp = load([mic_proc_dir fname]);

close all;
figure(1); set(gcf,'pos',[10 40 800 900],'color','w')
set(gca,'position',[.05 .05 .9 .9],'units','normalized')
hold on;
pb_mics=plot3(bpp.mic_loc(:,1),bpp.mic_loc(:,2),bpp.mic_loc(:,3),...
  '.k','markerfacecolor','k');
if plot_mic_nums
  text(bpp.mic_loc(:,1),bpp.mic_loc(:,2),bpp.mic_loc(:,3),...
    num2str((1:length(bpp.mic_loc))'),'color','r')
end

view(3); axis equal; grid on;