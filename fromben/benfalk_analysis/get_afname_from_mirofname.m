%get audio_fname from miro_fname
%returns 0 if not found
function afname=get_afname_from_mirofname(mirofname)

adir='..\mic_data\';
[~,~,raw] = xlsread('small_space_beampattern_analysis_log.xlsx',6);

C=strsplit(mirofname,{'_','.'});
cam=['c' C{2}];
trlindx=str2double( C{3} );

coli=find(strcmp(raw(1,:),cam));
rowi=find(cellfun(@(c) isequal(c,trlindx), raw(:,coli) ) );

sprow=find(~cellfun(@sum,cellfun(@isnan, raw(1:rowi,3),'uniformoutput',0)),1,'last');
daterow=find(~cellfun(@sum,cellfun(@isnan, raw(1:rowi,1),'uniformoutput',0)),1,'last');
bandrow=find(~cellfun(@sum,cellfun(@isnan, raw(1:rowi,4),'uniformoutput',0)),1,'last');

species=lower(raw{sprow,3});
date=datestr(raw(daterow,1),'yyyymmdd');
band=num2str(raw{bandrow,4});
trlnum=num2str(raw{rowi,11},'%.2i');

afname=strjoin({species,date,band,trlnum,'mic_data.mat'},'_');
if ~exist([adir afname],'file')
  afname=0;
end