function multifreq_data = get_multifreq_data(path,processed_files,freq_wanted)

for iB=1:length(processed_files)
    bat_proc_file = processed_files(iB).name;
    fprintf('File: %s\n',bat_proc_file);
    
    data = load(fullfile(path.base_path,path.data_path,bat_proc_file));
    if exist(fullfile(path.base_path,path.data_path,[bat_proc_file(1:end-4)...
        '_checked.mat']),'file')
      bpp_checked=load(fullfile(path.base_path,path.data_path,[bat_proc_file(1:end-4)...
        '_checked.mat']));
      data.proc.chk_good_call = bpp_checked.proc.chk_good_call;
      data.proc.ch_ex=bpp_checked.proc.ch_ex;
    end
    good_call_idx = find(data.proc.chk_good_call);
    

    for iC = good_call_idx'
        iC_save = find(iC==good_call_idx);
        fprintf('Call %02d\n',iC);
        
        for iF=1:length(freq_wanted)
            [call_dB,az,el,ch_include_idx] = get_call_azel_dB_data(data,freq_wanted(iF),iC);  % get good mic index
%             [vq,vq_norm,azq,elq] = interp_bp(az(ch_include_idx)/180*pi,el(ch_include_idx)/180*pi,call_dB(ch_include_idx),'rbf');  % use the first frequency data for finding ellipse center
            multifreq_data.call_dB{iB,iF}(iC_save,:) = call_dB;
            multifreq_data.call_dB_norm{iB,iF}(iC_save,:) = call_dB - max(call_dB);
            multifreq_data.az{iB,iF}(iC_save,:) = az;
            multifreq_data.el{iB,iF}(iC_save,:) = el;
            multifreq_data.ch_include_idx{iB,iF}(iC_save,:) = ch_include_idx;
%             multifreq_data.vq{iB,iF}{iC_save} = vq;
%             multifreq_data.vq_norm{iB,iF}{iC_save} = vq_norm;
%             multifreq_data.azq{iB,iF}{iC_save} = azq;
%             multifreq_data.elq{iB,iF}{iC_save} = elq;
        end
    end
end
