% 2015 11 08  Separate raw mic signal from detect file

DATA_DIR = 'E:\0_WJLEE_DATA\Beampattern Rousettus\20141219_rousettus\20141219_rousettus_mic_matfile';
fname = dir(fullfile(DATA_DIR,'*_detect.mat'));

for iF=1:length(fname)
    disp(fname(iF).name);
    m = matfile(fullfile(DATA_DIR,fname(iF).name),'Writable',true);
    m.sig = [];
    m.fs = [];
end