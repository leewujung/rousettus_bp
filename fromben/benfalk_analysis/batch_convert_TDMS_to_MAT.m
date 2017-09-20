% 2014 09 26  Batch convert TDMS data into MAT file

% sec_rec=625002/250e3;
% fs = 250e3;
% fs = 1e6;

GEN_DIR = 'E:\rousettus_closeup\';

folder = dir(GEN_DIR);
for iFF=3:length(folder)
  DATA_DIR = [GEN_DIR,filesep,folder(iFF).name];
  SAVE_DIR = [DATA_DIR,'_matfile'];
  files = dir([DATA_DIR,'/*_nodaqmx.tdms']);
  if ~exist(SAVE_DIR,'dir')
    mkdir(SAVE_DIR);
  end
  
  for iF=1:length(files)
    [ConvertedData,ConvertVer,ChanNames,GroupNames,ci]=convertTDMS(0,[DATA_DIR,'/',files(iF).name]);
    sig = [ConvertedData.Data.MeasuredData(:).Data];
    sig = sig-repmat(mean(sig,1),size(sig,1),1);
    xx= strcmp({ConvertedData.Data.MeasuredData(3).Property(:).Name},'wf_increment');
    fs = round(1/ConvertedData.Data.MeasuredData(3).Property(xx).Value);
    xx=find(strcmp({ConvertedData.Data.MeasuredData(3).Property(:).Name},'wf_start_time'));
    start_time=ConvertedData.Data.MeasuredData(3).Property(xx).Value;
    S = strsplit(files(iF).name,'_nodaqmx');
    save([SAVE_DIR,'/',S{1},'_mic_data.mat'],'sig','fs','start_time');
    
%     delete([DATA_DIR,'/',files(iF).name]);
    fff = strsplit([DATA_DIR,'/',files(iF).name],'.');
    fff = fff{1};
%     delete([fff,'.tdms_index']);
  end
end

