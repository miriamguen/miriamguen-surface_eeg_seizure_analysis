function score = get_local_score(EEG, onset, offset)
% onset and offset times will be in secondes 
EEG_local = EEG;
srate = EEG.srate;
EEG_local.data = EEG.data(:,floor(onset*srate):floor(offset*srate));

EEG_local.pnts =  length(EEG_local.data);
EEG_local.times = 1000*([0:EEG_local.pnts]/srate); 
EEG_local = pop_iclabel(EEG_local, 'default');
score = score_comp(EEG_local);
end