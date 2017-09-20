% 2014 09 26  Batch convert TDMS data into MAT file
% 2015 11 02  Grab fs from the tdms file

addpath C:\Users\wlee76\Dropbox\0_CODE\MATLAB\humphreysb-ConvertTDMS-878b34e

GEN_DIR = 'F:\Anand_Rousettus\Cindydata';

folder = dir(GEN_DIR);
for iFF=3:length(folder)
    DATA_DIR = [GEN_DIR,filesep,folder(iFF).name];
    SAVE_DIR = [DATA_DIR,'_matfile'];
    files = dir([DATA_DIR,'/*_nodaqmx.tdms']);
    if ~exist(SAVE_DIR,'dir')
        mkdir(SAVE_DIR);
    end
    
    for iF=1:length(files)
        [ConvertedData,ConvertVer,ChanNames,GroupNames,ci]=convertTDMS(0,[DATA_DIR,filesep,files(iF).name]);
        sig = [ConvertedData.Data.MeasuredData(:).Data];
        fs = 1/ConvertedData.Data.MeasuredData(end).Property(3).Value;  % tdms file records time increment
        sig = sig-repmat(mean(sig,1),size(sig,1),1);
        S = strsplit(files(iF).name,'_nodaqmx');
        save([SAVE_DIR,'/',S{1},'_mic_data.mat'],'sig','fs');
        
        delete([DATA_DIR,'/',files(iF).name]);
        fff = strsplit([DATA_DIR,'/',files(iF).name],'.');
        fff = fff{1};
        delete([fff,'.tdms_index']);
    end
end

