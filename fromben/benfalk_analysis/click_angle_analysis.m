%analysis of click angles with distance and with release time

clear;

species_date='rousettus_20150825';

%loads:
%  angle_diff,bat_dist2mics,beam_dirs_az,beam_dirs_el,callnums,voc_t,bat_pos,trials,beam_I
load(['extracted_beam_dirs_' species_date '.mat'])

PI_all=[]
for tt=1:length(angle_diff)  
  PI=diff(voc_t{tt})';
  pairs = PI < .03;
  
  pair_indx=sort([find(pairs) find(pairs)+1]);
  bd2m_pair{tt,1}=mean([bat_dist2mics{tt}(pairs) bat_dist2mics{tt}(find(pairs)+1)],2);
  
  voc_t_pairs{tt,1}=mean([voc_t{tt}(pairs) voc_t{tt}(find(pairs)+1)],2);
  
  bat_pos_pairs{tt,1}=bat_pos{tt}(pairs,:);
  
  good_pos_pairs{tt,1}=bat_pos_pairs{tt}(:,1)>.45;
  
  PI_all = [PI_all; PI'];
end

figure(30); clf;
t=0:.01:.25;
N=histcounts(PI_all,t);
plot(t(1:end-1) +.005/2,N)

AD=cell2mat(angle_diff')';
GPP=cell2mat(good_pos_pairs);


figure(1); clf; 
histogram(AD()./pi*180,10)
xlabel('Angle between pairs of clicks (deg)')


fnames={trials.name};
spl_fnames=cellfun(@(c) strsplit(c,'_'),fnames,'uniformoutput',0);
bands = cellfun(@(c) c{3},spl_fnames,'uniformoutput',0);
figure(4); clf; set(gcf,'color','w')
figure(7); clf; set(gcf,'color','w')
bats=unique(bands);
for bb=1:length(bats)
  bat_tt=strcmp(bands,bats{bb});
  ADbat=cell2mat(angle_diff(bat_tt)')';
  GPPbat=cell2mat(good_pos_pairs(bat_tt));
  Ibat=cell2mat(beam_I(bat_tt));
  
  figure(4);
  subaxis(length(bats),1,bb,'m',.05,'mb',.075,'sv',.05,'ml',.1)
  histogram(ADbat(:)./pi*180,10)
  ylabel(['Bat: ' bats{bb}]);
  if bb==length(bats)
    xlabel('Angle between pairs of clicks (deg)')
  end
  bat_angle_diff(bb) = mean(ADbat(:)/pi*180);
  
  figure(7);
  subaxis(length(bats),1,bb,'m',.05,'mb',.075,'sv',.05,'ml',.1)
  histogram(Ibat(:),10)
  ylabel(['Bat: ' bats{bb}]);
  if bb==length(bats)
    xlabel('Click Intensity')
  end
  Ibat_mean(bb)=mean(Ibat);
end
mean(bat_angle_diff)
std(bat_angle_diff)./sqrt(length(bats))



BD2MP = cell2mat(bd2m_pair);
figure(2), clf;
scatter(BD2MP(:),AD(:)/pi*180,'+')
xlabel('Distance btw. bat and nearest mic (m)')
ylabel('Click angle (deg)')
[b,bint,r,rint,stats]=regress(AD,[ones(size(AD)) BD2MP])
%R2 statistic, the F statistic and its p value, and an estimate of the error variance.


BPP=cell2mat(bat_pos_pairs);
figure(3), clf;
bubsizes = round( (min(AD(:)) :10/180*pi: max(AD))./pi*180 );
legentry=cell(size(bubsizes));
hold on
for ind = 1:numel(bubsizes)
   bubleg(ind) = plot(0,0,'ko','markersize',sqrt(bubsizes(ind)),'MarkerFaceColor','k');
   set(bubleg(ind),'visible','off')
   legentry{ind} = num2str(bubsizes(ind));
end
h=scatter3(BPP(:,1),BPP(:,2),BPP(:,3),AD(:)./pi*180,'k','filled');
alpha(h,.5)
legend(legentry,'location','southoutside','orientation','horizontal')
axis equal; view(2), grid on;


BD2M = cell2mat(bat_dist2mics);
BI=cell2mat(beam_I);
figure(5), clf; set(gcf,'color','w')
scatter(BD2M,BI,'+')
xlabel('Distance btw. bat and nearest mic (m)')
ylabel('Click Intensity')
LinearModel.fit(BD2M,BI)

figure(6), clf; set(gcf,'color','w')
histogram(BI)
xlabel('Click Intensity')




% figure(20),clf;
% VTP=cell2mat(voc_t_pairs);
% scatter(VTP(:),AD(:),'+')
% [b,bint,r,rint,stats]=regress(AD,[ones(size(AD)) VTP])

