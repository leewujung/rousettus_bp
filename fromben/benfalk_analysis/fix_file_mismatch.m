clear

comb_path='..\mic_data\';
combfiles=dir([comb_path 'eptesicus_20150824_LB62_*mic_data.mat']);
combfnames={combfiles.name}';

cal_path='..\mic_recordings\20150824_calib_mic\eptesicus_LB62_matfile\';
calfiles=dir([cal_path 'eptesicus_LB62*_mic_data.mat']);
calfnames={calfiles.name}';

pet_path='..\mic_recordings\20150824\eptesicus_LB62_matfile\';
petfiles=dir([pet_path 'eptescus_LB62*_mic_data.mat']);
petfnames={petfiles.name}';

[~,~,a]=xlsread('small_space_beampattern_analysis_log.xlsx','Cam_mic_all');

mic_nums=a(:,11);
cal_nums=a(:,12);

indx=239:247;

for k=1:length(indx)
  ffcal=find(strcmp(calfnames,...
    ['eptesicus_LB62' num2str(cal_nums{indx(k)}) '_mic_data.mat']));
  ffpet=find(strcmp(petfnames,...
    ['eptescus_LB62_' num2str(mic_nums{indx(k)}) '_mic_data.mat']));
  
  cal_data=load([cal_path calfnames{ffcal}]);
  y=resample(cal_data.sig,cal_data.fs/4,cal_data.fs);
  y=[zeros(1,size(y,2)); y];
  
  
  pet_data=load([pet_path petfnames{ffpet} ]);
  
  if length(y) ~= length(pet_data.sig)
    y = cal_data.sig;
  end
  
  if length(y) == length(pet_data.sig)
    
    pet_data.sig(:,end+1:end+size(y,2))=y;

    ffcomb = strcmp(combfnames,...
      ['eptesicus_20150824_LB62_' num2str(mic_nums{indx(k)},'%2.2d') '_mic_data.mat']);
    
    sig=pet_data.sig;
    copyfile([comb_path combfnames{ffcomb}],...
      [comb_path 'fixed\'])
    save([comb_path 'fixed\' combfnames{ffcomb}],'sig','-append')
  end
end


















clear

%create file list where cal_indx doesn't match mic_indx

[~,~,a]=xlsread('small_space_beampattern_analysis_log.xlsx','Cam_mic_all');

mic_nums=a(:,11);
cal_nums=a(:,12);

mic_nan=cellfun(@isnan,mic_nums(2:end),'uniformoutput',0);
cal_exist=~cellfun(@isnan,cal_nums(2:end));

unmatched=[];
for k=find(cal_exist)'
  if ~isequal(mic_nums{k+1},cal_nums{k+1})
    unmatched(end+1)=k+1;
  end
  
end




clear

cal12fn='..\mic_recordings\20150824_calib_mic\eptesicus_LB62_matfile\eptesicus_LB6212_mic_data.mat';
cal12=load(cal12fn);
% pet12fn='..\mic_recordings\20150824\eptesicus_LB62_matfile\eptescus_LB62_12_mic_data.mat';
% pet12=load(pet12fn);
comb12_fn='..\mic_data\eptesicus_20150824_LB62_12_mic_data.mat';
comb12=load(comb12_fn);

size(cal12.sig)
size(comb12.sig)

figure(1);
subplot(2,1,1)
plot((1:length(cal12.sig))/cal12.fs,cal12.sig)
axis tight
subplot(2,1,2)
plot((1:length(comb12.sig))/comb12.fs,comb12.sig(:,33))
axis tight