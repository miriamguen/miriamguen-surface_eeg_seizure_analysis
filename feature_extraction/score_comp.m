function scores = score_comp(EEG)
% EEG - EEGLAB struct after runing mara and iclabel on it
feature_params
% create componant relaiability score - R - what is that?
ic_label = EEG.etc.ic_classification.ICLabel.classifications;
try
    mara_artifact_prob = EEG.etc.ic_classification.MARA.info.posterior_artefactprob;
    mara = (1- mara_artifact_prob)';
catch ME
    % warning(['this file is to short for mara, using brain label insted ' ME.identifier])
    mara = ic_label(:,1) + ic_label(:,7);
end
brain    = ic_label(:,1);   %(ic_label(:,1)>ic_brain_th)*2;
muscle   = ic_label(:,2);  %ic_label(:,2)<ic_th;
eye      = ic_label(:,3);  %ic_label(:,3)<ic_th;
heart    = ic_label(:,4);  %ic_label(:,4)<ic_th;
line     = ic_label(:,5);   %ic_label(:,5)<ic_th;
chan     = ic_label(:,6);  %ic_label(:,6)<ic_th;
other    = ic_label(:,7);  %ic_label(:,7)<ic_th;
scores   = table(mara, brain,muscle,eye, heart, line, chan,other);
% compute the componant score, the highest the score the mor likelly it is
% to be brain
scores.score = mean(scores{:,1:2},2);
end
