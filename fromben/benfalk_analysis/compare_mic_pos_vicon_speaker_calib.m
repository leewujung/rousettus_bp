err=nan(length(mic_infer_loc),1);
indx=nan(length(mic_infer_loc),1);
for k=1:length(mic_infer_loc)
  D=distance(mic_infer_loc(k,:),mic_pos);
  [err(k),indx(k)]=min(D);
end


SPM=cellfun(@mean,speaker_pos,'uniformoutput',0);

SPMlks=cell2mat(SPM)

text(mic_infer_loc(:,1),mic_infer_loc(:,2),mic_infer_loc(:,3),...
  num2str([1:length(mic_infer_loc)]'))