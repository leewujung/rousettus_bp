%get vicon from audio_fname
%returns 0 if not found
function [fname,f_exist]=get_vfname_from_afname(afname)
[~,~,raw] = xlsread('small_space_beampattern_analysis_log.xlsx',6);

C=strsplit(afname,{'_','.'});
date_name=[C{2}(6) '/' C{2}(7:8) '/' C{2}(1:4)];
species_name=C{1};
band=C{3};
trl_num=str2double( C{4} );

rowi_date=find(cellfun(@(c) isequal(c,date_name), raw(:,1) ) );
if length(rowi_date)>1
  next_rowi_date = find( cellfun(@(c) isequal(ischar(c),1), ...
    raw(rowi_date(end)+1:end,1) ) ,1)+rowi_date(end)-1;
  rowi_date=rowi_date(1);
else
  next_rowi_date = find(cellfun(@(c) isequal(ischar(c),1), raw(rowi_date+1:end,1) ) ,1)+rowi_date-1;
end
rowi_species=find( cellfun(@(c) isequal(lower(c),species_name), ...
  raw(rowi_date:next_rowi_date,3) ) ,1)+rowi_date-1;
rowi_band=find(cellfun(@(c) ~isempty(find(strfind(c,band), 1)) || ...
  isequal(c,str2double(band)), ...
  raw(rowi_species:next_rowi_date,4) ) )+rowi_species-1;
rowi_trial=find(cellfun(@(c) isequal(c,trl_num), ...
  raw(rowi_band:next_rowi_date,11) ) )+rowi_band-1;
rowi_trial=rowi_trial(1);

mdir='..\bat_pos\';

fname=[strjoin({species_name,datestr(date_name,'yyyymmdd'),band,...
  num2str(raw{rowi_trial,10},'%2.2i'),'bat_pos.mat'},'_')];
f_exist=exist([mdir fname],'file');